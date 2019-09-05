//
// Created by pkua on 05.09.2019.
//

#include "PovrayMode.h"
#include "../PackingGenerator.h"

PovrayMode::PovrayMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("<file in>"));

    this->packingFilename = positionalArguments[0];
}

void PovrayMode::run() {
    Packing packing;
    packing.restore(this->packingFilename);
    PackingGenerator::toPovray(packing, this->params.surfaceSize, nullptr, this->packingFilename + ".pov");
}
