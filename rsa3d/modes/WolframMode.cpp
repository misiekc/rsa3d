//
// Created by pkua on 05.09.2019.
//

#include "WolframMode.h"
#include "../utils/Assertions.h"
#include "../PackingGenerator.h"

WolframMode::WolframMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() < 1 || positionalArguments.size() > 3)
        die(arguments.formatUsage("<file in> (use periodic image = false) (bc expand fraction = 0.1)"));

    this->packingFilename = positionalArguments[0];

    if (positionalArguments.size() == 1 || positionalArguments[1] == "false")
        this->isPeriodicImage = false;
    else if (positionalArguments[1] == "true")
        this->isPeriodicImage = true;
    else
        die("(use periodic image) must be empty, 'true' or 'false'");

    if (positionalArguments.size() == 3) {
        this->bcExpandFraction = std::stod(positionalArguments[2]);
        ValidateMsg(this->bcExpandFraction >= 0 && this->bcExpandFraction < 0.5,
                    "BC expand fraction must be in [0, 0.5) range");
    } else {
        this->bcExpandFraction = 0.1;
    }
}

void WolframMode::run() {
    Packing packing;
    packing.restore(this->packingFilename);
    PackingGenerator::toWolfram(packing, this->params.surfaceSize, nullptr, this->isPeriodicImage,
                                this->bcExpandFraction, this->packingFilename + ".nb");
}
