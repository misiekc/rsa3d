//
// Created by pkua on 05.09.2019.
//

#include "ExclusionZonesMode.h"
#include "../analizator/ExclusionZoneVisualizer.h"

void ExclusionZonesMode::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("[packing file]"));

    this->packingFilename = positionalArguments[0];
}

void ExclusionZonesMode::run() {
    ExclusionZoneVisualizer::main(this->params, this->packingFilename);
}

void ExclusionZonesMode::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << arguments.formatUsage("[packing file]") << std::endl;
    out << std::endl;
    out << "Currently it is broken, so I am not sure what it does." << std::endl;
}
