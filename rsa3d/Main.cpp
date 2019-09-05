#include "PackingGenerator.h"
#include "analizator/Analyzer.h"
#include "analizator/ExclusionZoneVisualizer.h"
#include "shape/ShapeFactory.h"
#include "utils/Assertions.h"
#include "ProgramArguments.h"
#include "Simulation.h"
#include "simulation_modes/DefaultSimulation.h"
#include "simulation_modes/BoundariesSimulation.h"
#include "simulation_modes/AccuracySimulation.h"
#include "simulation_modes/DensitySimulation.h"


#include <sys/wait.h>
#include <iomanip>
#include <fstream>
#include <memory>



void makeDatFileForPackingsInDirectory(const Parameters *params, const std::string &dirName) {
    std::string sFile = params->getPackingSignature() + ".dat";
    std::ofstream dataFile(sFile);
    if (!dataFile)
        die("Cannot open file " + sFile + " to store packing info");
    dataFile.precision(std::numeric_limits<double>::digits10 + 1);

    auto filenames = PackingGenerator::findPackingsInDir(dirName);
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

/*int simulate(const ProgramArguments &arguments) {
    std::string sFile = params.getPackingSignature() + ".dat";
    std::ofstream file(sFile);
    if (!file)
        die("Cannot open file " + sFile + " to store packing info");
    file.precision(std::numeric_limits<double>::digits10 + 1);

    std::size_t seed = params.from;
    if (params.generatorProcesses > 1) {
        while (seed < params.from + params.collectors) {
            std::size_t i;
            for (i = 0; ((i < params.generatorProcesses) &&
                         ((seed + i) < (params.from + params.collectors))); i++) {
                int pid = fork();
                if (pid < 0) {
                    std::cout << "fork problem" << std::endl;
                    i--;
                    continue;
                }
                if (pid == 0) {
                    runSingleSimulation(static_cast<int>(seed + i), &params, file);
                    return 1;
                }
                if (pid > 0) {
                    continue;
                }
            }
            for (std::size_t j = 0; j < i; j++) {
                wait(nullptr);
                seed++;
            }
        }
    } else {
        for (std::size_t i = 0; i < params.collectors; i++) {
            runSingleSimulation(static_cast<int>(params.from + i), &params, file);
        }
    }
    file.close();
    return 1;
}*/


void test(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 2)
        die(arguments.formatUsage("<file in> <max time>"));

    Packing packing;
    packing.restore(positionalArguments[0]);
    PackingGenerator pg(1, &arguments.getParameters());
    pg.testPacking(packing, std::stod(positionalArguments[1]));
}

void debug(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("<packing generator file>"));

    Parameters params = arguments.getParameters();
    const char *cfile = positionalArguments[0].c_str();
    std::string filename(cfile);
    char buf[20];
    sprintf(buf, "%.0f", pow(params.surfaceSize, params.surfaceDimension));
    std::string size(buf);

    std::ifstream file(filename, std::ios::binary);
    if (!file)
        die("Cannot open file " + filename + " to restore packing generator");

    PackingGenerator pg(0, &params);
    pg.restore(file);
    pg.run();
}

void dat(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("<directory>"));

    std::string directory = positionalArguments[0];
    makeDatFileForPackingsInDirectory(&arguments.getParameters(), directory);
}

void analyze(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() < 1 || positionalArguments.size() > 2)
        die(arguments.formatUsage("<directory> (correlations range = 10; 0 - no corr output)"));

    double corrRange = (positionalArguments.size() >= 2) ? std::stod(positionalArguments[1]) : 10.0;
    Validate(corrRange >= 0);
    Analyzer an(&arguments.getParameters());
    an.analyzePackingsInDirectory(positionalArguments[0], 0.01, 1.0, corrRange);
}

void povray(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("<file in>"));

    std::string file(positionalArguments[0]);
    Packing packing;
    packing.restore(file);
    PackingGenerator::toPovray(packing, arguments.getParameters().surfaceSize, nullptr, file + ".pov");
}

void wolfram(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() < 1 || positionalArguments.size() > 3)
        die(arguments.formatUsage("<file in> (use periodic image = false) (bc expand fraction = 0.1)"));

    bool isPeriodicImage{};
    if (positionalArguments.size() == 1 || positionalArguments[1] == "false")
        isPeriodicImage = false;
    else if (positionalArguments[1] == "true")
        isPeriodicImage = true;
    else
        die("(use periodic image) must be empty, 'true' or 'false'");

    double bcExpandFraction{};
    if (positionalArguments.size() == 3) {
        bcExpandFraction = std::stod(positionalArguments[2]);
        ValidateMsg(bcExpandFraction >= 0 && bcExpandFraction < 0.5, "BC expand fraction must be in [0, 0.5) range");
    } else {
        bcExpandFraction = 0.1;
    }

    std::string file(positionalArguments[0]);
    Packing packing;
    packing.restore(file);
    PackingGenerator::toWolfram(packing, arguments.getParameters().surfaceSize, nullptr, isPeriodicImage,
                                bcExpandFraction, file + ".nb");
}

void bc_expand(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() < 1 || positionalArguments.size() > 2)
        die(arguments.formatUsage("<file in> (file out = file in)"));

    std::string fileIn(positionalArguments[0]);
    std::string fileOut = (positionalArguments.size() == 1 ? fileIn : positionalArguments[1]);
    Packing packing;
    packing.restore(fileIn);
    packing.expandOnPBC(arguments.getParameters().surfaceSize, 0.1);
    packing.store(fileOut);
}

void exclusion_zones(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 2)
        die(arguments.formatUsage("<packing file> <output file>"));

    std::string packingFile(positionalArguments[0]);
    std::string outputFile(positionalArguments[1]);
    ExclusionZoneVisualizer::main(arguments.getParameters(), packingFile, outputFile);
}

int main(int argc, char **argv) {
    std::unique_ptr<ProgramArguments> arguments;
    try {
        arguments = std::make_unique<ProgramArguments>(argc, argv);
    } catch (InvalidArgumentsException &e) {
        die(e.what());
    }

    Parameters params = arguments->getParameters();
    ShapeFactory::initShapeClass(params.particleType, params.particleAttributes);

    std::unique_ptr<Simulation> simulationMode;
    std::string mode = arguments->getMode();

    // Simulation modes - they all generate a number of packings, operate on them, create *.dat file and optionally
    // save additional information
    if (mode == "simulate")
        simulationMode = std::make_unique<DefaultSimulation>(*arguments);
    else if (mode == "boundaries")
        simulationMode = std::make_unique<BoundariesSimulation>(*arguments);
    else if (mode == "accuracy")
        simulationMode = std::make_unique<AccuracySimulation>(*arguments);
    else if (mode == "density")
        simulationMode = std::make_unique<DensitySimulation>(*arguments);

    // Other modes
    else if (mode == "test")
        test(*arguments);
    else if (mode == "debug")
        debug(*arguments);
    else if (mode == "dat")
        dat(*arguments);
    else if (mode == "analyze")
        analyze(*arguments);
    else if (mode == "povray")
        povray(*arguments);
    else if (mode == "wolfram")
        wolfram(*arguments);
    else if (mode == "bc_expand")
        bc_expand(*arguments);
    else if (mode == "exclusion_zones")
        exclusion_zones(*arguments);
    else
        die("Unknown mode: " + mode);

    if (simulationMode)
        simulationMode->run();

    return EXIT_SUCCESS;
}
