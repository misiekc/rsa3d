//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_BOUNDARIESSIMULATION_H
#define RSA3D_BOUNDARIESSIMULATION_H


#include "Simulation.h"
#include "../ProgramArguments.h"

class BoundariesSimulation : public Simulation {
private:
    std::size_t targetNumberOfParticles{};
    std::size_t particleCounter{};

public:
    void initializeForArguments(const ProgramArguments &arguments) override;
    void printHelp(std::ostream &out, const ProgramArguments &arguments) override;

    bool continuePackingGeneration(std::size_t packingIndex, const Packing &packing) override;
    void postProcessPacking(Packing &packing) override;

    void postProcessSimulation() override { }
};


#endif //RSA3D_BOUNDARIESSIMULATION_H
