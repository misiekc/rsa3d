//
// Created by pkua on 05.09.2019.
//

#include <fstream>

#include "DatFileGenerationMode.h"
#include "../PackingGenerator.h"

DatFileGenerationMode::DatFileGenerationMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("<directory>"));

    this->dirName = positionalArguments[0];
}

void DatFileGenerationMode::run() {
    std::string sFile = this->params.getPackingSignature() + ".dat";
    std::ofstream dataFile(sFile);
    if (!dataFile)
        die("Cannot open file " + sFile + " to store packing info");
    dataFile.precision(std::numeric_limits<double>::digits10 + 1);

    auto filenames = PackingGenerator::findPackingsInDir(this->dirName);
    for (const auto &filename : filenames) {
        int no1 = lastIndexOf(filename, '_');
        int no2 = lastIndexOf(filename, '.');
        std::string seed = filename.substr(no1 + 1, no2 - no1 - 1);

        Packing packing;
        packing.restore(filename);

        dataFile << seed << "\t" << packing.back()->no << "\t" << packing.back()->time << std::endl;
        std::cout << ".";
        std::cout.flush();
    }
    std::cout << std::endl;
}
