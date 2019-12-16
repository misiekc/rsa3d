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

void PackingTestMode::printHelp(std::ostream &out, const std::string &cmd) {
    out << "Usage: " << cmd << " test [packing file] [max time]" << std::endl;
    out << std::endl;
    out << "It loads the packing from [packing file] and tries to add a new particle up to [max time]" << std::endl;
    out << "dimensionless time. It can be used to check if the packing is saturated. Note, that it" << std::endl;
    out << "requires to be invoked with the same parameters as were used to generate the packing." << std::endl;
}
