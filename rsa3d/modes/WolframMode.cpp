//
// Created by pkua on 05.09.2019.
//

#include "WolframMode.h"
#include "../utils/Assertions.h"
#include "../PackingGenerator.h"

void WolframMode::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() < 1 || positionalArguments.size() > 2)
        die(arguments.formatUsage("<file in> (use periodic image = false)"));

    this->packingFilename = positionalArguments[0];

    if (positionalArguments.size() == 1 || positionalArguments[1] == "false")
        this->isPeriodicImage = false;
    else if (positionalArguments[1] == "true")
        this->isPeriodicImage = true;
    else
        die("(use periodic image) must be empty, 'true' or 'false'");
}

void WolframMode::run() {
    Packing packing;
    packing.restore(this->packingFilename);
    PackingGenerator::toWolfram(packing, this->params.surfaceSize, nullptr, this->isPeriodicImage,
                                this->packingFilename + ".nb");
}
