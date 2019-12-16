//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_SIMULATION_H
#define RSA3D_SIMULATION_H


#include "../Packing.h"
#include "../Parameters.h"
#include "../ProgramMode.h"

class Simulation : public ProgramMode {
private:
    Packing runSingleSimulation(int seed, std::ofstream &dataFile);

protected:
    virtual bool continuePackingGeneration(std::size_t packingIndex, const Packing &packing) = 0;
    virtual void postProcessPacking(Packing &packing) = 0;
    virtual void postProcessSimulation() = 0;

public:
    void printHelp(std::ostream &out, const ProgramArguments &arguments) override;

    void run() final;
};


#endif //RSA3D_SIMULATION_H
