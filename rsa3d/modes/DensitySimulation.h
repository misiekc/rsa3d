//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_DENSITYSIMULATION_H
#define RSA3D_DENSITYSIMULATION_H


#include "Simulation.h"
#include "../ProgramArguments.h"

class DensitySimulation : public Simulation {
private:
    double targetPackingFraction{};

public:
    explicit DensitySimulation(const ProgramArguments &arguments);

    bool continuePackingGeneration(std::size_t packingIndex, const Packing &packing) override;
    void postProcessPacking(Packing &packing) override;

    void postProcessSimulation() override { }
};


#endif //RSA3D_DENSITYSIMULATION_H
