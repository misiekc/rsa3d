//
// Created by pkua on 05.09.2019.
//

#include "PackingTestMode.h"
#include "../utils/Assertions.h"
#include "../PackingGenerator.h"

void PackingTestMode::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 2)
        die(arguments.formatUsage("<file in> <max time>"));

    this->packingFilename = positionalArguments[0];
    this->maxTime = std::stod(positionalArguments[1]);
    Validate(this->maxTime > 0);
}

void PackingTestMode::run() {
    Packing packing;
    packing.restore(this->packingFilename);
    PackingGenerator pg(1, &this->params);
    pg.testPacking(packing, this->maxTime);
}
