//
// Created by pkua on 05.09.2019.
//

#include "BoundariesSimulation.h"

void BoundariesSimulation::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("[number of particles]"));

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

void BoundariesSimulation::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << arguments.formatUsage("[number of particles]") << std::endl;
    out << std::endl;
    Simulation::printHelp(out, arguments);
    out << "It perform as many simulation as needed to obtain a total count of particles greater than" << std::endl;
    out << "[number of particles]." << std::endl;
}
