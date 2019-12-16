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

void DefaultSimulation::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << arguments.formatUsage("") << std::endl;
    out << std::endl;
    Simulation::printHelp(out, arguments);
    out << "It perform as many simulation as specified by 'collectors' parameter and does not perform any" << std::endl;
    out << "additional operations." << std::endl;
}
