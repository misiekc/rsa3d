#include <sys/wait.h>
#include <memory>

#include "shape/ShapeFactory.h"
#include "ProgramArguments.h"
#include "modes/DefaultSimulation.h"
#include "modes/BoundariesSimulation.h"
#include "modes/AccuracySimulation.h"
#include "modes/DensitySimulation.h"
#include "modes/PackingTestMode.h"
#include "modes/DebugMode.h"
#include "modes/DatFileGenerationMode.h"
#include "modes/AnalyzeMode.h"
#include "modes/PovrayMode.h"
#include "modes/BCExpandMode.h"
#include "modes/ExclusionZonesMode.h"
#include "modes/WolframMode.h"


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
