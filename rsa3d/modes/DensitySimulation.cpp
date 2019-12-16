//
// Created by pkua on 05.09.2019.
//

#include "DensitySimulation.h"

void DensitySimulation::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (!positionalArguments.empty())
        die(arguments.formatUsage(""));
}

bool DensitySimulation::continuePackingGeneration(std::size_t packingIndex, const Packing &packing) {
    return packingIndex < this->params.collectors - 1;
}

void DensitySimulation::postProcessPacking(Packing &packing) {
    double currentPackingFraction = packing.getParticlesVolume(this->params.surfaceDimension)
                                    / this->params.sufraceVolume();

    if (this->params.maxDensity > currentPackingFraction) {
        std::cout << "[DensitySimulation::postProcessPacking] Target packing fraction: ";
        std::cout << this->params.maxDensity << " > generated packing fraction: ";
        std::cout << currentPackingFraction << ". Keeping original packing." << std::endl;
    } else {
        std::cout << "[DensitySimulation::postProcessPacking] Initial packing fraction      : ";
        std::cout << currentPackingFraction << ", reducing... " << std::flush;
        double targetVolume = this->params.maxDensity * this->params.sufraceVolume();
        packing.reducePackingVolume(targetVolume, this->params.surfaceDimension);
        std::cout << "done." << std::endl;

        double actualPackingFraction = packing.getParticlesVolume(this->params.surfaceDimension)
                                       / this->params.sufraceVolume();
        std::cout << "[DensitySimulation::postProcessPacking] Target packing fraction       : ";
        std::cout << this->params.maxDensity << std::endl;
        std::cout << "[DensitySimulation::postProcessPacking] Actual packing fraction       : ";
        std::cout << actualPackingFraction << std::endl;
        std::cout << "[DensitySimulation::postProcessPacking] Last particle adsorption time : ";
        std::cout << packing.back()->time << std::endl;
    }
}
