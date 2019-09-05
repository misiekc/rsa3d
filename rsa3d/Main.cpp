#include "PackingGenerator.h"
#include "analizator/Analyzer.h"
#include "analizator/ExclusionZoneVisualizer.h"
#include "shape/ShapeFactory.h"
#include "utils/Assertions.h"
#include "ProgramArguments.h"
#include "modes/Simulation.h"
#include "modes/DefaultSimulation.h"
#include "modes/BoundariesSimulation.h"
#include "modes/AccuracySimulation.h"
#include "modes/DensitySimulation.h"


#include <sys/wait.h>
#include <iomanip>
#include <fstream>
#include <memory>

namespace {
    class PackingTestMode : public ProgramMode {
    private:
        std::string packingFilename;
        double maxTime{};

    public:
        explicit PackingTestMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
            std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
            if (positionalArguments.size() != 2)
                die(arguments.formatUsage("<file in> <max time>"));

            this->packingFilename = positionalArguments[0];
            this->maxTime = std::stod(positionalArguments[1]);
            Validate(this->maxTime > 0);
        }

        void run() override {
            Packing packing;
            packing.restore(this->packingFilename);
            PackingGenerator pg(1, &this->params);
            pg.testPacking(packing, this->maxTime);
        }
    };


    class DebugMode : public ProgramMode {
    private:
        std::string packingGeneratorFilename;

    public:
        explicit DebugMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
            std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
            if (positionalArguments.size() != 1)
                die(arguments.formatUsage("<packing generator file>"));

            this->packingGeneratorFilename = positionalArguments[0];
        }

        void run() override {
            std::ifstream file(this->packingGeneratorFilename, std::ios::binary);
            if (!file)
                die("Cannot open file " + this->packingGeneratorFilename + " to restore packing generator");

            PackingGenerator pg(0, &params);
            pg.restore(file);
            pg.run();
        }
    };

    class DatFileGenerationMode : public ProgramMode {
    private:
        std::string dirName;

    public:
        explicit DatFileGenerationMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
            std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
            if (positionalArguments.size() != 1)
                die(arguments.formatUsage("<directory>"));

            this->dirName = positionalArguments[0];
        }

        void run() override {
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
    };

    class AnalyzeMode : public ProgramMode {
    private:
        std::string dirName;
        double corrRange{};
    public:
        explicit AnalyzeMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
            std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
            if (positionalArguments.size() < 1 || positionalArguments.size() > 2)
                die(arguments.formatUsage("<directory> (correlations range = 10; 0 - no corr output)"));

            this->dirName = positionalArguments[0];
            this->corrRange = (positionalArguments.size() >= 2) ? std::stod(positionalArguments[1]) : 10.0;
            Validate(this->corrRange >= 0);
        }

        void run() override {
            Analyzer an(&this->params);
            an.analyzePackingsInDirectory(this->dirName, 0.01, 1.0, this->corrRange);
        }
    };

    class PovrayMode : public ProgramMode {
    private:
        std::string packingFilename;

    public:
        explicit PovrayMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
            std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
            if (positionalArguments.size() != 1)
                die(arguments.formatUsage("<file in>"));

            this->packingFilename = positionalArguments[0];
        }

        void run() override {
            Packing packing;
            packing.restore(this->packingFilename);
            PackingGenerator::toPovray(packing, this->params.surfaceSize, nullptr, this->packingFilename + ".pov");
        }
    };

    class WolframMode : public ProgramMode {
    private:
        std::string packingFilename;
        bool isPeriodicImage{};
        double bcExpandFraction{};

    public:
        explicit WolframMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
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

        void run() override {
            Packing packing;
            packing.restore(this->packingFilename);
            PackingGenerator::toWolfram(packing, this->params.surfaceSize, nullptr, this->isPeriodicImage,
                                        this->bcExpandFraction, this->packingFilename + ".nb");
        }
    };

    class BCExpandMode : public ProgramMode {
    private:
        std::string packingInFilename;
        std::string packingOutFilename;

    public:
        explicit BCExpandMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
            std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
            if (positionalArguments.size() < 1 || positionalArguments.size() > 2)
                die(arguments.formatUsage("<file in> (file out = file in)"));

            this->packingInFilename = positionalArguments[0];
            this->packingOutFilename = (positionalArguments.size() == 1 ?
                                        this->packingInFilename : positionalArguments[1]);
        }

        void run() override {
            Packing packing;
            packing.restore(this->packingInFilename);
            packing.expandOnPBC(this->params.surfaceSize, 0.1);
            packing.store(this->packingOutFilename);
        }
    };

    class ExclusionZonesMode : public ProgramMode {
    private:
        std::string packingFilename;
        std::string outputFilename;

    public:
        explicit ExclusionZonesMode(const ProgramArguments &arguments) : ProgramMode(arguments.getParameters()) {
            std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
            if (positionalArguments.size() != 2)
                die(arguments.formatUsage("<packing file> <output file>"));

            this->packingFilename = positionalArguments[0];
            this->outputFilename = positionalArguments[1];
        }

        void run() override {
            ExclusionZoneVisualizer::main(this->params, this->packingFilename, this->outputFilename);
        }
    };
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

int main(int argc, char **argv) {
    std::unique_ptr<ProgramArguments> arguments;
    try {
        arguments = std::make_unique<ProgramArguments>(argc, argv);
    } catch (InvalidArgumentsException &e) {
        die(e.what());
    }

    Parameters params = arguments->getParameters();
    ShapeFactory::initShapeClass(params.particleType, params.particleAttributes);

    std::unique_ptr<ProgramMode> programMode;
    std::string mode = arguments->getMode();

    // Simulation modes - they all generate a number of packings, operate on them, create *.dat file and optionally
    // save additional information
    if (mode == "simulate")
        programMode = std::make_unique<DefaultSimulation>(*arguments);
    else if (mode == "boundaries")
        programMode = std::make_unique<BoundariesSimulation>(*arguments);
    else if (mode == "accuracy")
        programMode = std::make_unique<AccuracySimulation>(*arguments);
    else if (mode == "density")
        programMode = std::make_unique<DensitySimulation>(*arguments);

    // Other modes
    else if (mode == "test")
        programMode = std::make_unique<PackingTestMode>(*arguments);
    else if (mode == "debug")
        programMode = std::make_unique<DebugMode>(*arguments);
    else if (mode == "dat")
        programMode = std::make_unique<DatFileGenerationMode>(*arguments);
    else if (mode == "analyze")
        programMode = std::make_unique<AnalyzeMode>(*arguments);
    else if (mode == "povray")
        programMode = std::make_unique<PovrayMode>(*arguments);
    else if (mode == "wolfram")
        programMode = std::make_unique<WolframMode>(*arguments);
    else if (mode == "bc_expand")
        programMode = std::make_unique<BCExpandMode>(*arguments);
    else if (mode == "exclusion_zones")
        programMode = std::make_unique<ExclusionZonesMode>(*arguments);
    else
        die("Unknown mode: " + mode);

    programMode->run();
    return EXIT_SUCCESS;
}
