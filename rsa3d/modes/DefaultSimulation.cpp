//
// Created by pkua on 05.09.2019.
//

#include "DefaultSimulation.h"

DefaultSimulation::DefaultSimulation(const ProgramArguments &arguments) : Simulation(arguments.getParameters()) {
    if (!arguments.getPositionalArguments().empty())
        die(arguments.formatUsage(""));
}

bool DefaultSimulation::continuePackingGeneration(std::size_t packingIndex, const Packing &packing) {
    return packingIndex < this->params.collectors - 1;
}
