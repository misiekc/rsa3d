//
// Created by ciesla on 11/14/22.
//

#include <cstdlib>
#include <fstream>
#include "DomainAnalyzer.h"
#include "../shape/shapes/DiscreteOrientationsShape2_1.h"
#include "../PackingGenerator.h"
#include "../shape/shapes/Rectangle.h"
#include "../../statistics/LogPlot.h"

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
        std::cout << ".";
    }

    std::cout << std::endl;
    std::cout << dirName << "\t" << sq / counter << "\t"
              << sqrt((sq2 / counter - sq * sq / (counter * counter)) / counter) << std::endl;
}

void DomainAnalyzer::exportPacking(const std::string &fileIn, const std::string &fileOut){
    Packing packing;
    packing.restore(fileIn);
    std::ofstream file(fileOut);
    for(const RSAShape *s : packing){
        file << s->getPosition()[0] << "\t" << s->getPosition()[1] << "\t";
        if (s->getOrientation()[0]==0){
            file << Rectangle::longer << "\t" << Rectangle::shorter;
        }else{
            file << Rectangle::shorter << "\t" << Rectangle::longer;
        }
        file << std::endl;
    }
}

void DomainAnalyzer::dividePackingIntoDomains(const Packing &packing, std::vector<Domain> &domains){
    bool nextLoop;
    std::vector<const RSAShape *> copy = packing.getVector();
    while (copy.size()>0) {
        Domain d;
        d.addShape(copy[0]);
        copy.erase(copy.begin());
        nextLoop = true;
        while (nextLoop) {
            nextLoop = false;
            for (size_t i = 0; i < copy.size(); i++) {
                if (d.testShape(copy[i])) {
                    d.addShape(copy[i]);
                    copy.erase(copy.begin() + i);
                    i--;
                    nextLoop = true;
                    break;
                }
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
    for (const auto &packingPath: packingPaths) {
        Packing packing;
        packing.restore(packingPath);
        std::vector<Domain> domains;
        DomainAnalyzer::dividePackingIntoDomains(packing, domains);
        for (Domain d: domains){
            domainSizes.push_back(d.size());
        }
        std::cout << ".";
    }
    std::cout << std::endl;
    std::sort(domainSizes.begin(), domainSizes.end());
    LogPlot pDomains(domainSizes[0], domainSizes[domainSizes.size()-1], 100);
    for(size_t size: domainSizes) {
        pDomains.add(size);
    }
    printHistogram(pDomains, dirName + "_domains.txt");
}
