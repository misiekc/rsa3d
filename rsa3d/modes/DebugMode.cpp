//
// Created by pkua on 05.09.2019.
//

#include <fstream>
#include "DebugMode.h"
#include "../PackingGenerator.h"

void DebugMode::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("<packing generator file>"));

    this->packingGeneratorFilename = positionalArguments[0];
}

void DebugMode::run() {
    std::ifstream file(this->packingGeneratorFilename, std::ios::binary);
    if (!file)
        die("Cannot open file " + this->packingGeneratorFilename + " to restore packing generator");

    PackingGenerator pg(0, &params);
    pg.restore(file);
    pg.run();
}
