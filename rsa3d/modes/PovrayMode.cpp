//
// Created by pkua on 05.09.2019.
//

#include "PovrayMode.h"
#include "../PackingGenerator.h"

void PovrayMode::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() < 1 || positionalArguments.size() > 2)
        die(arguments.formatUsage("[file in] (use periodic image = false)"));

    this->packingFilename = positionalArguments[0];

    if (positionalArguments.size() == 1 || positionalArguments[1] == "false")
        this->isPeriodicImage = false;
    else if (positionalArguments[1] == "true")
        this->isPeriodicImage = true;
    else
        die("(use periodic image) must be empty, 'true' or 'false'");
}

void PovrayMode::run() {
    Packing packing;
    packing.restore(this->packingFilename);
    PackingGenerator::toPovray(packing, this->params.surfaceSize, nullptr, this->isPeriodicImage,
                               this->packingFilename + ".pov");
}

void PovrayMode::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << arguments.formatUsage("(use periodic image = false)") << std::endl;
    out << std::endl;
    out << "It draws a packing using povray. If [use periodic image] is true, the packing will be" << std::endl;
    out << "expanded according to the boundary conditions so that the resulting image can be placed" << std::endl;
    out << "repeatedly on a grid and the borders will match, Note, that it requires to be invoked with" << std::endl;
    out << "the same parameters as were used to generate the packing." << std::endl;
}
