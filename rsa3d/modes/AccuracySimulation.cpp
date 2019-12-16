//
// Created by pkua on 05.09.2019.
//

#include "AccuracySimulation.h"

void AccuracySimulation::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 2)
        die(arguments.formatUsage("[accuracy] [out file]"));
    this->targetAccuracy = std::stod(positionalArguments[0]);
    Validate(this->targetAccuracy > 0);
    this->outputFilename = positionalArguments[1];
}

bool AccuracySimulation::continuePackingGeneration(std::size_t packingIndex, const Packing &packing) {
    return packingIndex < 2 || this->meanPackingFraction.error > this->targetAccuracy;
}

void AccuracySimulation::postProcessPacking(Packing &packing) {
    this->packingFractions.push_back(static_cast<double>(packing.size()) / this->params.sufraceVolume());
    this->meanPackingFraction.calculateFromSamples(this->packingFractions);
    std::cout << "[AccuracySimulation::postProcessPacking] Current accuracy : ";
    std::cout << this->meanPackingFraction.error << std::endl;
    std::cout << "[AccuracySimulation::postProcessPacking] Target accuracy  : ";
    std::cout << this->targetAccuracy << std::endl;
}

void AccuracySimulation::postProcessSimulation() {
    std::ofstream file(this->outputFilename, std::ios_base::app);
    if (!file)
        die("Cannot open " + this->outputFilename + " file to store packing fraction");
    file.precision(std::numeric_limits<double>::digits10 + 1);
    file << this->meanPackingFraction.value << "\t" << this->meanPackingFraction.error << "\t";
    file << this->params.particleType << "\t" << this->params.particleAttributes << std::endl;
}

void AccuracySimulation::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << arguments.formatUsage("[accuracy] [file out]") << std::endl;
    out << std::endl;
    Simulation::printHelp(out, arguments);
    out << "It perform as many simulation as needed for packing fraction mean error to be smaller than" << std::endl;
    out << "[accuracy]. Moreover, it appends to [file out] a line containing particle attributes and" << std::endl;
    out << "packing fraction with error." << std::endl;
}
