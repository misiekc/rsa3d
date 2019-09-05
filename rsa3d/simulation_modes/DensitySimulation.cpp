//
// Created by pkua on 05.09.2019.
//

#include "DensitySimulation.h"

DensitySimulation::DensitySimulation(const ProgramArguments &arguments) : Simulation(arguments.getParameters()) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("<target packing fraction>"));
    this->targetPackingFraction = std::stod(positionalArguments[0]);
    Validate(this->targetPackingFraction > 0.0 && this->targetPackingFraction <= 1.0);
}

bool DensitySimulation::continuePackingGeneration(std::size_t packingIndex, const Packing &packing) {
    return packingIndex < this->params.collectors - 1;
}

void DensitySimulation::postProcessPacking(Packing &packing) {
    double currentPackingFraction = packing.getParticlesVolume(this->params.surfaceDimension)
                                    / this->params.sufraceVolume();

    if (this->targetPackingFraction > currentPackingFraction) {
        std::cout << "[DensitySimulation::postProcessPacking] Target packing fraction: ";
        std::cout << this->targetPackingFraction << " > generated packing fraction: ";
        std::cout << currentPackingFraction << ". Keeping original packing." << std::endl;
    } else {
        std::cout << "[DensitySimulation::postProcessPacking] Initial packing fraction      : ";
        std::cout << currentPackingFraction << ", reducing... " << std::flush;
        double targetVolume = this->targetPackingFraction * this->params.sufraceVolume();
        packing.reducePackingVolume(targetVolume, this->params.surfaceDimension);
        std::cout << "done." << std::endl;

        double actualPackingFraction = packing.getParticlesVolume(this->params.surfaceDimension)
                                       / this->params.sufraceVolume();
        std::cout << "[DensitySimulation::postProcessPacking] Target packing fraction       : ";
        std::cout << this->targetPackingFraction << std::endl;
        std::cout << "[DensitySimulation::postProcessPacking] Actual packing fraction       : ";
        std::cout << actualPackingFraction << std::endl;
        std::cout << "[DensitySimulation::postProcessPacking] Last particle adsorption time : ";
        std::cout << packing.back()->time << std::endl;
    }
}
