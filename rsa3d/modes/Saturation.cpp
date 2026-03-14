//
// Created by pkua on 05.09.2019.
//

#include <chrono>
#include <fstream>
#include <unistd.h>
#include <filesystem>

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
    std::filesystem::path path(this->packingFilename);
    std::error_code ec;
    std::vector<std::string> filenames{};
    if (std::filesystem::is_directory(path, ec) ){
        auto allFilenames = PackingGenerator::findPackingsInDir(this->packingFilename);
        for (auto sfile : allFilenames) {
            if (sfile.find(".ns.")!=std::string::npos)
                filenames.push_back(sfile);
        }
    }else {
        filenames.push_back(this->packingFilename);
    }
    for (auto sfile : filenames) {
        packing.restore(sfile);
        std::size_t collector = this->getCollectorNumber(sfile);
        unsigned long seedOrigin = this->getSeedOrigin();
        unsigned long seed = seedOrigin + collector;
        PackingGenerator pg(seed, collector, &this->params);
        bool bSaturated = pg.run(&packing);
        if (params.storePackings) {
            std::string spath = sfile.substr(0, sfile.rfind(std::filesystem::path::preferred_separator)+1);
            std::string sfilename = spath + pg.getPackingFilename(bSaturated);
            packing.store(sfilename);
        }
        std::cout << std::endl << std::flush;
    }
}

unsigned long Saturation::getCollectorNumber(std::string sfile) {
    size_t i0 = sfile.rfind("_")+1;
    size_t i1 = sfile.find(".ns.bin");
    std::string sindex = sfile.substr(i0, i1-i0);
    return static_cast<unsigned long>(atol(sindex.c_str()));
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

