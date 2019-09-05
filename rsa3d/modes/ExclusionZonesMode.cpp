//
// Created by pkua on 05.09.2019.
//

#include "ExclusionZonesMode.h"
#include "../analizator/ExclusionZoneVisualizer.h"

ExclusionZonesMode::ExclusionZonesMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 2)
        die(arguments.formatUsage("<packing file> <output file>"));

    this->packingFilename = positionalArguments[0];
    this->outputFilename = positionalArguments[1];
}

void ExclusionZonesMode::run() {
    ExclusionZoneVisualizer::main(this->params, this->packingFilename, this->outputFilename);
}
