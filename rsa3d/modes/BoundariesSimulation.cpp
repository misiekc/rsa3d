//
// Created by pkua on 05.09.2019.
//

#include "BoundariesSimulation.h"

BoundariesSimulation::BoundariesSimulation(const ProgramArguments &arguments) : Simulation(arguments.getParameters()) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("<number of particles>"));

    this->targetNumberOfParticles = std::stoul(positionalArguments[0]);
    Validate(this->targetNumberOfParticles > 0);
}

bool BoundariesSimulation::continuePackingGeneration(std::size_t packingIndex, const Packing &packing) {
    return this->particleCounter < this->targetNumberOfParticles;
}

void BoundariesSimulation::postProcessPacking(Packing &packing) {
    this->particleCounter += packing.size();
    std::cout << "[BoundariesSimulation::postProcessPacking] Fraction of particles generated: ";
    std::cout << static_cast<double>(this->particleCounter) / this->targetNumberOfParticles << std::endl;
}
