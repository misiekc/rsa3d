//
// Created by pkua on 05.09.2019.
//

#include "BCExpandMode.h"
#include "../ProgramArguments.h"
#include "../Packing.h"

void BCExpandMode::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() < 1 || positionalArguments.size() > 2)
        die(arguments.formatUsage("<file in> (file out = file in)"));

    this->packingInFilename = positionalArguments[0];
    this->packingOutFilename = (positionalArguments.size() == 1 ? this->packingInFilename : positionalArguments[1]);
}

void BCExpandMode::run() {
    Packing packing;
    packing.restore(this->packingInFilename);
    packing.expandOnPBC(this->params.surfaceSize);
    packing.store(this->packingOutFilename);
}
