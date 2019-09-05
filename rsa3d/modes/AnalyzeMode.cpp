//
// Created by pkua on 05.09.2019.
//

#include "AnalyzeMode.h"
#include "../utils/Assertions.h"
#include "../analizator/Analyzer.h"

AnalyzeMode::AnalyzeMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() < 1 || positionalArguments.size() > 2)
        die(arguments.formatUsage("<directory> (correlations range = 10; 0 - no corr output)"));

    this->dirName = positionalArguments[0];
    this->corrRange = (positionalArguments.size() >= 2) ? std::stod(positionalArguments[1]) : 10.0;
    Validate(this->corrRange >= 0);
}

void AnalyzeMode::run() {
    Analyzer an(&this->params);
    an.analyzePackingsInDirectory(this->dirName, 0.01, 1.0, this->corrRange);
}
