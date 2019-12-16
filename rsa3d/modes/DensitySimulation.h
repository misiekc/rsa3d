//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_DENSITYSIMULATION_H
#define RSA3D_DENSITYSIMULATION_H


#include "Simulation.h"
#include "../ProgramArguments.h"

class DensitySimulation : public Simulation {
public:
    void initializeForArguments(const ProgramArguments &arguments) override;

    void printHelp(std::ostream &out, const std::string &cmd) override;

    bool continuePackingGeneration(std::size_t packingIndex, const Packing &packing) override;
    void postProcessPacking(Packing &packing) override;

    void postProcessSimulation() override { }
};


#endif //RSA3D_DENSITYSIMULATION_H
