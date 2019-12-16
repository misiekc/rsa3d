//
// Created by pkua on 05.09.2019.
//

#include "AnalyzeMode.h"
#include "../utils/Assertions.h"
#include "../analizator/Analyzer.h"

void AnalyzeMode::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
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

void AnalyzeMode::printHelp(std::ostream &out, const std::string &cmd) {
    out << "Usage: " << cmd << " analyze [directory] (correlacions range = 10; 0 - no corr output)" << std::endl;
    out << std::endl;
    out << "It searches [directory] for '*.bin' files and analyzes them. It produces a few files, each" << std::endl;
    out << "with prefix '[directory]_': kinetics plot (suffix 'nvt'), ASF plot ('asf'), order parameters," << std::endl;
    out << "if available for the shape used ('order'), and prints on the standard output parameters such" << std::endl;
    out << "as packing fraction and d, all labeled. Note, that it requires to be invoked with the same" << std::endl;
    out << "parameters as were used to generate the packings." << std::endl;
}
