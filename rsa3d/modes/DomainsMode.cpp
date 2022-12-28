//
// Created by Michal Ciesla on 16.11.2022.
//

#include "DomainsMode.h"
#include "../utils/Assertions.h"
#include "../analizator/DomainAnalyzer.h"
#include <filesystem>

void DomainsMode::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() < 1)
        die(arguments.formatUsage("name - directory name with packings to analyze or filename to draw"));

    this->dirFile = positionalArguments[0];
}

void DomainsMode::run() {
    ValidateMsg(this->params.particleType.rfind("DiscreteOrientations", 0)==0, "Only DiscreteOrientations shape are supported by this mode");
    ValidateMsg(this->params.particleAttributes.rfind("2 ", 0)==0, "Only two discrete orientations are supported by this mode");
    DomainAnalyzer analyzer(this->params);
    if (std::filesystem::is_directory((this->dirFile))) {
        analyzer.analyzeOrderDirectory(this->dirFile);
        analyzer.analyzeDomains(this->dirFile);
    }else if(std::filesystem::is_regular_file(this->dirFile)) {
        Packing packing;
        packing.restore(this->dirFile);
        analyzer.toWolfram(packing, this->dirFile + ".nb");
    }
}

void DomainsMode::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << arguments.formatUsage("dirFile [fileOut]") << std::endl;
    out << std::endl;
    out << "It searches [directory] for '*.bin' files and analyzes domains. It produces a few file" << std::endl;
    out << "with prefix '[directory]_domains.txt' and prints on the standard output order parameter." << std::endl;
}
