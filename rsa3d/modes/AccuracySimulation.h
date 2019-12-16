//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_ACCURACYSIMULATION_H
#define RSA3D_ACCURACYSIMULATION_H


#include <fstream>

#include "Simulation.h"
#include "../ProgramArguments.h"
#include "../utils/Quantity.h"

class AccuracySimulation : public Simulation {
private:
    double targetAccuracy{};
    std::string outputFilename;
    std::vector<double> packingFractions;
    Quantity meanPackingFraction;

public:
    void initializeForArguments(const ProgramArguments &arguments) override;
    bool continuePackingGeneration(std::size_t packingIndex, const Packing &packing) override;
    void postProcessPacking(Packing &packing) override;
    void postProcessSimulation() override;
};


#endif //RSA3D_ACCURACYSIMULATION_H
