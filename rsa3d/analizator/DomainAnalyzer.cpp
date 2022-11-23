//
// Created by Michal Ciesla on 11/14/22.
//

#include <cstdlib>
#include <fstream>
#include <unistd.h>
#include "DomainAnalyzer.h"
#include "../shape/shapes/DiscreteOrientationsShape2_1.h"
#include "../PackingGenerator.h"

Parameters DomainAnalyzer::params;

void DomainAnalyzer::init(const Parameters &p){
    DomainAnalyzer::params = p;
}

double DomainAnalyzer::analyzeOrder(const Packing &packing) {
    size_t nh=0, nv=0;
    for (const RSAShape *si :packing) {
        if (si->getOrientation()[0] == 0)
            nh++;
        else
            nv++;
    }
    return 1 - 2.0 * static_cast<double>(std::min(nv, nh) / static_cast<double>(nv + nh));
}

void DomainAnalyzer::analyzeOrderDirectory(const std::string &dirName) {
    auto packingPaths = PackingGenerator::findPackingsInDir(dirName);
    double sq = 0.0, sq2 = 0.0;
    int counter = 0;
    std::cout << "[DomainAnalyzer::analyzePackingsInDirectory] " << std::flush;
    for (const auto &packingPath: packingPaths) {
        Packing packing;
        packing.restore(packingPath);
        double q = DomainAnalyzer::analyzeOrder(packing);
        sq += q;
        sq2 += q * q;
        counter++;
        std::cout << "." << std::flush;
    }
    std::cout << std::endl;
    std::cout << dirName << "\t" << sq / counter << "\t"
              << sqrt((sq2 / counter - sq * sq / (counter * counter)) / counter) << std::endl;
}

void DomainAnalyzer::addNeighboursToQueue(std::vector<const RSAShape*> &neighbours, std::unordered_set<const RSAShape *> &queue){
    for (const RSAShape *s : neighbours){
        if (queue.find(s)==queue.end())
            queue.insert(s);
    }
}

void DomainAnalyzer::dividePackingIntoDomains(const Packing &packing, std::vector<const Domain *> &domains){
    std::unordered_set<const RSAShape *> copy;
    copy.insert(packing.getVector().begin(), packing.getVector().end());
    RSAPeriodicBC bc(DomainAnalyzer::params.surfaceSize);

    while (copy.size()>0) {
        Domain* d = new Domain(DomainAnalyzer::params);
        std::unordered_set<const RSAShape *> analyzeStack;
        std::vector<const RSAShape *> neighbours;
        NeighbourGrid<const RSAShape> ng(DomainAnalyzer::params.surfaceDimension, DomainAnalyzer::params.surfaceSize, RSAShape::getNeighbourListCellSize());
        for(const RSAShape *s: copy) {
            ng.add(s, s->getPosition());
        }
        const RSAShape *s = *(copy.begin());
        copy.erase(s);
        d->addShape(s);
        ng.remove(s, s->getPosition());
        ng.getNeighbours(&neighbours, s->getPosition());
        DomainAnalyzer::addNeighboursToQueue(neighbours, analyzeStack);

        while (analyzeStack.size()>0) {
            const RSAShape *s = *analyzeStack.begin();
            analyzeStack.erase(s);
            if (d->testShape(s)) {
                d->addShape(s);
                std::unordered_set<const RSAShape *>::iterator it = copy.find(s);
                if (it!=copy.end()) {
                    copy.erase(it);
                    ng.remove(s, s->getPosition());
                }else {
                    std::cout << "Something wrong" << std::endl;
                }
                ng.getNeighbours(&neighbours, s->getPosition());
                DomainAnalyzer::addNeighboursToQueue(neighbours, analyzeStack);
                if (d->size()%1000 !=0)
                    continue;
            }
        }
        domains.push_back(d);
    }
}

void DomainAnalyzer::printHistogram(Plot& plot, const std::string &filename){
    double **plotPoints = new double*[plot.size()];
    for(int i=0; i<plot.size(); i++){
        plotPoints[i] = new double[2];
    }
    plot.getAsHistogramPoints(plotPoints);

    std::ofstream file(filename);

    for (int i = 0; i < plot.size()-1; i++) {
        if (plotPoints[i][1] != 0) {
            file << plotPoints[i][0] << "\t" << plotPoints[i][1];
            file << std::endl;
        }
    }
    file.close();

    for(int i=0; i<plot.size(); i++){
        delete[] plotPoints[i];
    }
    delete[] plotPoints;
}


void DomainAnalyzer::analyzeDomains(const std::string &dirName) {
    std::vector<size_t> domainSizes;
    auto packingPaths = PackingGenerator::findPackingsInDir(dirName);
    std::cout << "[DomainAnalyzer::analyzeDomains] " << std::flush;
    std::string rawFilename = dirName + "_domains.raw";
    if (access(rawFilename.c_str(), F_OK)==0) { // if exists file the analysis is obsolete
        std::ifstream rawFile(rawFilename);
        size_t size;
        while (!rawFile.eof()) {
            rawFile >> size;
            _OMP_CRITICAL(domainSizes)
            domainSizes.push_back(size);
        }
        rawFile.close();
    }else{ // performing domain analysis
        _OMP_PARALLEL_FOR
        for (const auto &packingPath: packingPaths) {
            Packing packing;
            packing.restore(packingPath);
            std::vector<const Domain *> domains;
            DomainAnalyzer::dividePackingIntoDomains(packing, domains);
            std::ofstream rawFile(rawFilename);
            for (const Domain *d: domains){
                size_t size = d->size();
                delete d;
                rawFile << size << std::endl;
                _OMP_CRITICAL(domainSizes)
                domainSizes.push_back(size);
            }
            rawFile.close();
        }
        std::cout << "." << std::flush;
    }
    std::cout << std::endl;

    Plot pDomains(domainSizes[0], domainSizes[domainSizes.size()-1], 50);
    for(size_t size: domainSizes) {
        pDomains.add(size);
    }
    printHistogram(pDomains, dirName + "_domains.txt");
}

void DomainAnalyzer::toWolfram(Packing packing, const std::string &filename) {
    std::vector<std::string> colors = {"Red", "Green", "Blue", "Black", "Gray", "Cyan", "Magenta", "Yellow", "Brown", "Orange", "Pink", "Purple"};
    size_t colorsSize = colors.size();
    std::vector<const Domain *> domains;
    DomainAnalyzer::dividePackingIntoDomains(packing, domains);
    size_t i = 0;
    std::ofstream file(filename);
    file << "Graphics[{" << colors[i%colorsSize];
    for (const Domain* d: domains){
        i++;
        for (const RSAShape *s : d->shapes)
            file << ", " << std::endl << s->toWolfram();
        file << ", " << colors[i%colorsSize];
    }
    file << "}";
    file << "]" << std::endl;

    file.close();
}

