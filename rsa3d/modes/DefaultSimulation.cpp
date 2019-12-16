//
// Created by pkua on 05.09.2019.
//

#include "DefaultSimulation.h"

void DefaultSimulation::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    if (!arguments.getPositionalArguments().empty())
        die(arguments.formatUsage(""));
}

bool DefaultSimulation::continuePackingGeneration(std::size_t packingIndex, const Packing &packing) {
    return packingIndex < this->params.collectors - 1;
}

void DefaultSimulation::printHelp(std::ostream &out, const std::string &cmd) {
    out << "Usage: " << cmd << " simulate" << std::endl;
    out << std::endl;
    out << "Performs a standard simulation for a given parameters without anything additional." << std::endl;
}
