//
// Created by pkua on 05.09.2019.
//

#include <chrono>
#include <fstream>
#include <unistd.h>

#include "Saturation.h"
#include "../PackingGenerator.h"


void Saturation::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("[file in]"));

    this->packingFilename = positionalArguments[0];
}

void Saturation::run() {
    Packing packing;
    packing.restore(this->packingFilename);
    std::size_t collector = params.from;
    unsigned long seedOrigin = this->getSeedOrigin();
    unsigned long seed = seedOrigin + collector;
    PackingGenerator pg(seed, collector, &this->params);
    pg.run(&packing);
    if (params.storePackings)
        packing.store(pg.getPackingFilename());

}

unsigned long Saturation::getSeedOrigin() const {
    if (this->params.seedOrigin == "random")
        return std::random_device{}();
    else
        return std::stoul(this->params.seedOrigin);
}

void Saturation::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << arguments.formatUsage("[file in]") << std::endl;
    out << std::endl;
    out << "It loads a given packing and tries to saturate it." << std::endl;
}

