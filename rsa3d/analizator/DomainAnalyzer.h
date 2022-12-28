//
// Created by Michal Ciesla on 11/14/22.
//

#ifndef RSA3D_DOMAINANALYZER_H
#define RSA3D_DOMAINANALYZER_H

#include <unordered_set>
#include "../Packing.h"
#include "../shape/shapes/Rectangle.h"
#include "../shape/shapes/DiscreteOrientationsShape2_1.h"
#include "../../statistics/Plot.h"
#include "../Parameters.h"
#include "../NeighbourGrid.h"
#include "../boundary_conditions/PeriodicBC.h"
#include "../../statistics/LogPlot.h"
#include "../boundary_conditions/FreeBC.h"
#include "Domain.h"

class DomainAnalyzer {
private:
    Parameters params;
    double analyzeOrder(const Packing &packing);
    void addNeighboursToQueue(std::vector<const RSAShape*> &neighbours, std::unordered_set<const RSAShape *> &queue);

    template<typename InputIterator>
    void getComponents(InputIterator begin, InputIterator end, RSABoundaryConditions *bc, std::vector<const Domain *> &domains);

    bool isPercolating(const Domain &domain);
    void printHistogram(Plot &plot, const std::string &filename);
    void toWolfram(const std::vector<const Domain *> &domains, const std::string &filename);

public:
    DomainAnalyzer(const Parameters &params);
    void analyzeOrderDirectory(const std::string &dirName);
    void analyzeDomains(const std::string &dirName);
    void toWolfram(const Packing& packing, const std::string &filename);
};

#endif //RSA3D_DOMAINANALYZER_H
