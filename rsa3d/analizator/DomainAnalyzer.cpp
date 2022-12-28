//
// Created by Michal Ciesla on 11/14/22.
//

#include <fstream>
#include <unistd.h>
#include "DomainAnalyzer.h"
#include "../PackingGenerator.h"

DomainAnalyzer::DomainAnalyzer(const Parameters &p){
    this->params = p;
}

double DomainAnalyzer::analyzeOrder(const Packing &packing) {
    size_t nh=0, nv=0;
    for (const RSAShape *si :packing) {
        if (si->getOrientation()[0] == 0)
            nh++;
        else
            nv++;
    }
    return 1 - 2.0 * static_cast<double>(std::min(nv, nh)) / static_cast<double>(nv + nh);
}

void DomainAnalyzer::analyzeOrderDirectory(const std::string &dirName) {
    auto packingPaths = PackingGenerator::findPackingsInDir(dirName);
    double sq = 0.0, sq2 = 0.0;
    int counter = 0;
    std::cout << "[DomainAnalyzer::analyzeOrderDirectory] " << std::flush;
    for (const auto &packingPath: packingPaths) {
        Packing packing;
        packing.restore(packingPath);
        double q = this->analyzeOrder(packing);
        sq += q;
        sq2 += q * q;
        counter++;
        std::cout << "." << std::flush;
    }
    std::cout << std::endl;
    std::cout << dirName << "\t" << sq / counter << "\t"
              << sqrt((sq2 / counter - sq * sq / (counter * counter)) / counter) << std::endl;
}

void DomainAnalyzer::addNeighboursToQueue(std::vector<const RSAShape *> &neighbours,
                                          std::unordered_set<const RSAShape *> &queue) {
    for (const RSAShape *s : neighbours){
        if (queue.find(s)==queue.end())
            queue.insert(s);
    }
}
template<typename InputIterator>
void DomainAnalyzer::getComponents(InputIterator begin, InputIterator end, RSABoundaryConditions *bc, std::vector<const Domain *> &domains){
    std::unordered_set<const RSAShape *> copy;
    copy.insert(begin, end);

    while (!copy.empty()) {
        std::unordered_set<const RSAShape *> analyzeStack;
        std::vector<const RSAShape *> neighbours;
        NeighbourGrid<const RSAShape> ng(this->params.surfaceDimension, this->params.surfaceSize, RSAShape::getNeighbourListCellSize());
        for(const RSAShape *s: copy) {
            ng.add(s, s->getPosition());
        }
        const RSAShape *seed = *(copy.begin());
        copy.erase(seed);
        Domain *d = new Domain(this->params.surfaceSize, this->params.surfaceDimension, seed);
        ng.remove(seed, seed->getPosition());
        ng.getNeighbours(&neighbours, seed->getPosition());
        this->addNeighboursToQueue(neighbours, analyzeStack);

        while (!analyzeStack.empty()) {
            const RSAShape *s = *analyzeStack.begin();
            analyzeStack.erase(s);
            if (d->testAndAddShape(s, bc)) {
                auto it = copy.find(s);
                if (it!=copy.end()) {
                    copy.erase(it);
                    ng.remove(s, s->getPosition());
                }else {
                    std::cout << "Something wrong" << std::endl;
                }
                ng.getNeighbours(&neighbours, s->getPosition());
                this->addNeighboursToQueue(neighbours, analyzeStack);
                if (d->size()%1000 !=0)
                    continue;
            }
        }
        domains.push_back(d);
    }
}

bool DomainAnalyzer::isPercolating(const Domain &domain){
    RSAFreeBC fbc;
    std::vector<const Domain *> domains;
    this->getComponents(domain.getShapes().begin(), domain.getShapes().end(), &fbc, domains);
//    DomainAnalyzer::toWolfram(domains, "divided.nb");
    for (const Domain *d: domains) {
        double xMin = this->params.surfaceSize, yMin = this->params.surfaceSize, xMax = 0, yMax = 0;
        for (const RSAShape *s: d->getShapes()) {
            RSAVector pos = s->getPosition();
            if (pos[0] < xMin) xMin = pos[0];
            if (pos[1] < yMin) yMin = pos[1];
            if (pos[0] > xMax) xMax = pos[0];
            if (pos[1] > yMax) yMax = pos[1];
        }
        if (d->getOrientation() == 0) {
            if ((xMax - xMin) > this->params.surfaceSize - (Rectangle::longer + Rectangle::shorter))
                return true;
            if ((yMax - yMin) > this->params.surfaceSize - (Rectangle::shorter + Rectangle::shorter))
                return true;
        } else {
            if ((yMax - yMin) > this->params.surfaceSize - (Rectangle::longer + Rectangle::shorter))
                return true;
            if ((xMax - xMin) > this->params.surfaceSize - (Rectangle::shorter + Rectangle::shorter))
                return true;
        }
    }
    return false;
}

void DomainAnalyzer::printHistogram(Plot& plot, const std::string &filename){
    auto **plotPoints = new double*[plot.size()];
    for(size_t i=0; i<plot.size(); i++){
        plotPoints[i] = new double[2];
    }
    plot.getAsHistogramPoints(plotPoints);

    std::ofstream file(filename);

    for (size_t i = 0; i < plot.size()-1; i++) {
        if (plotPoints[i][1] != 0) {
            file << plotPoints[i][0] << "\t" << plotPoints[i][1];
            file << std::endl;
        }
    }
    file.close();

    for(size_t i=0; i<plot.size(); i++){
        delete[] plotPoints[i];
    }
    delete[] plotPoints;
}

void DomainAnalyzer::analyzeDomains(const std::string &dirName) {
    RSAPeriodicBC pbc(this->params.surfaceSize);
    std::vector<size_t> domainSizes;
    auto packingPaths = PackingGenerator::findPackingsInDir(dirName);
    std::cout << "[DomainAnalyzer::analyzeDomains] " << std::flush;
    std::string rawFilename = dirName + "_domains.raw";
    size_t percolating = 0;
    if (access(rawFilename.c_str(), F_OK)==0) { // if exists file the analysis is obsolete
        std::ifstream rawFile(rawFilename);
        size_t size;
        while (!rawFile.eof()) {
            rawFile >> size;
            domainSizes.push_back(size);
        }
        rawFile.close();
    }else{ // performing domain analysis
        _OMP_PARALLEL_FOR
        for (const auto &packingPath: packingPaths) {
            Packing packing;
            packing.restore(packingPath);
            std::vector<const Domain *> domains;
            this->getComponents(packing.getVector().begin(), packing.getVector().end(), &pbc, domains);
            for (const Domain *d: domains){
/*
                std::vector<const Domain*> v;
                v.push_back(d);
                this->toWolfram(v, packingPath + "_tested.nb");
*/
                size_t size = d->size();
                if (this->isPercolating(*d)) {
                    _OMP_CRITICAL(percolating)
                    percolating++;
                    this->toWolfram(domains, packingPath + "_percolation.nb");
                    break;
                }
                _OMP_CRITICAL(domainSizes)
                domainSizes.push_back(size);
            }
            for (const Domain *d: domains) {
                delete d;
            }
            std::cout << "." << std::flush;
        }
        std::sort(domainSizes.begin(), domainSizes.end());
        std::ofstream rawFile(rawFilename);
        for(size_t size : domainSizes)
            rawFile << size << std::endl;
        rawFile.close();
    }
    std::cout << std::endl;
    std::cout << static_cast<double>(percolating) / static_cast<double>(packingPaths.size());

    LogPlot pDomainsLog(static_cast<double>(domainSizes[0]), static_cast<double>(domainSizes[domainSizes.size()-1]), 50);
//    Plot pDomainsLin(static_cast<double>(domainSizes[0]), static_cast<double>(domainSizes[domainSizes.size()-1]), 50);
    for (size_t size: domainSizes) {
        pDomainsLog.add(static_cast<double>(size));
//        pDomainsLin.add(static_cast<double>(size));
    }

    printHistogram(pDomainsLog, dirName + "_domains.txt");
//    printHistogram(pDomainsLin, dirName + "_domains_lin.txt");
}

void DomainAnalyzer::toWolfram(const std::vector<const Domain *> &domains, const std::string &filename) {
    std::vector<std::string> colors = {"Red", "Green", "Blue", "Black", "Gray", "Cyan", "Magenta", "Yellow", "Brown", "Orange", "Pink", "Purple"};
    size_t colorsSize = colors.size();
    size_t i = 0;
    std::ofstream file(filename);
    file << "Graphics[{" << colors[i%colorsSize];
    for (const Domain* d: domains){
        i++;
        for (const RSAShape *s : d->getShapes())
            file << ", " << std::endl << s->toWolfram();
        file << ", " << colors[i%colorsSize];
    }
    file << "}";
    file << "]" << std::endl;

    file.close();
}

void DomainAnalyzer::toWolfram(const Packing& packing, const std::string &filename) {
    std::vector<const Domain *> domains;
    RSAPeriodicBC pbc(this->params.surfaceSize);
    this->getComponents(packing.getVector().begin(), packing.getVector().end(), &pbc, domains);
    this->toWolfram(domains, filename);
    for(const Domain *d : domains)
        delete d;
}
