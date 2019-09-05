//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_SIMULATION_H
#define RSA3D_SIMULATION_H


#include "Packing.h"
#include "Parameters.h"

class Simulation {
private:
    Packing runSingleSimulation(int seed, std::ofstream &dataFile);

protected:
    Parameters params;

    virtual bool continuePackingGeneration(std::size_t packingIndex, const Packing &packing) = 0;
    virtual void postProcessPacking(Packing &packing) = 0;
    virtual void postProcessSimulation() = 0;

public:
    explicit Simulation(Parameters params) : params{std::move(params)} { }

    void run();
};


#endif //RSA3D_SIMULATION_H
