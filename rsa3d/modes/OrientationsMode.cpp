//
// Created by pkua on 05.09.2019.
//

#include "OrientationsMode.h"
#include "../Packing.h"
#include "../analizator/OrientationAnalyzer.h"

void OrientationsMode::initializeForArguments(const ProgramArguments &arguments) {
    if (RSA_SPATIAL_DIMENSION!=3 || RSA_ANGULAR_DIMENSION!=3) {
        die("only for 3,3 shapes");
    }
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() < 1 || positionalArguments.size() > 2)
        die(arguments.formatUsage("[file in]"));
    this->packingFilename = positionalArguments[0];
}

void OrientationsMode::run() {
    Packing packing;
    packing.restore(this->packingFilename);
    OrientationAnalyzer::analyzeOrientations(packing);
}

void OrientationsMode::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << arguments.formatUsage("file in") << std::endl;
    out << std::endl;
    out << "TODO: orientations analyzer" << std::endl;
}
