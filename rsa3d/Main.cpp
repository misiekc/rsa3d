#include <memory>
#include <map>
#include <iterator>
#include <functional>

#include "shape/ShapeFactory.h"
#include "ProgramArguments.h"
#include "modes/DefaultSimulation.h"
#include "modes/Saturation.h"
#include "modes/BoundariesSimulation.h"
#include "modes/AccuracySimulation.h"
#include "modes/DensitySimulation.h"
#include "modes/PackingTestMode.h"
#include "modes/DebugMode.h"
#include "modes/DatFileGenerationMode.h"
#include "modes/AnalyzeMode.h"
#include "modes/PovrayMode.h"
#include "modes/BCExpandMode.h"
#include "modes/ContactFunctionMode.h"
#include "modes/ExclusionZonesMode.h"
#include "modes/WolframMode.h"
#include "modes/SurfaceAreaMode.h"
#include "modes/DomainsMode.h"


namespace {
    // This it the place where you want to add additional program modes. Just implement ProgramMode / Simulation
    // and add it here. The rest it already done for you ;)
    std::map<std::string, std::shared_ptr<ProgramMode>> programModes = {
            // Simulation modes - they all generate a number of packings, operate on them, create *.dat file and
            // optionally save additional information
            {"simulate", std::make_shared<DefaultSimulation>()},
            {"saturate", std::make_shared<Saturation>()},
            {"boundaries", std::make_shared<BoundariesSimulation>()},
            {"accuracy", std::make_shared<AccuracySimulation>()},
            {"density", std::make_shared<DensitySimulation>()},
            // Other modes
            {"test", std::make_shared<PackingTestMode>()},
            {"debug", std::make_shared<DebugMode>()},
            {"dat", std::make_shared<DatFileGenerationMode>()},
            {"analyze", std::make_shared<AnalyzeMode>()},
//            {"domains", std::make_shared<DomainsMode>()},
            {"povray", std::make_shared<PovrayMode>()},
            {"wolfram", std::make_shared<WolframMode>()},
            {"bc_expand", std::make_shared<BCExpandMode>()},
            {"contact_function", std::make_shared<ContactFunctionMode>()},
            {"exclusion_zones", std::make_shared<ExclusionZonesMode>()},
            {"surface_area", std::make_shared<SurfaceAreaMode>()}};

    void print_modes(std::ostream &out) {
        using modes_pair_t = decltype(programModes)::value_type;
        using namespace std::placeholders;
        std::transform(programModes.begin(), programModes.end(), std::ostream_iterator<std::string>(out, ", "),
                       std::bind(&modes_pair_t::first, _1));
    }

    void print_general_help(std::ostream &out, const ProgramArguments &arguments) {
        out << arguments.formatUsage("") << std::endl;
        out << std::endl;
        out << "At the beginning the program initializes parameters with default values. They can be changed using ";
        out << "flags:" << std::endl;
        out << "  -f [file] : reads the parameters from a given file begin a key=value config" << std::endl;
        out << "  -i : reads the parameters from a standard input, in the form as above" << std::endl;
        out << "  -[parameter name]=[parameter value] : anything else starting with - is treated as a parameter (no ";
        out << "spaces!)" << std::endl;
        out << "You can combine all the above methods of specifying parameters. The latter ones override the former.";
        out << std::endl << std::endl;
        out << "Available modes: ";
        print_modes(out);
        out << std::endl;
        out << "To see the help for a specific mode, use '" << arguments.getCmd() << " help [mode]'" << std::endl;
    }

    int handle_help_request(const ProgramArguments &arguments) {
        std::string mode = arguments.getMode();
        if (mode.empty()) {
            print_general_help(std::cout, arguments);
            return EXIT_SUCCESS;
        }

        if (programModes.find(mode) != programModes.end()) {
            programModes[mode]->printHelp(std::cout, arguments);
            return EXIT_SUCCESS;
        } else {
            std::cerr << "Unknown mode: " << mode << ". Type '" << arguments.getCmd() << " help' for help";
            std::cerr << std::endl;
            return EXIT_FAILURE;
        }
    }
}

int main(int argc, char **argv) {
    std::unique_ptr<ProgramArguments> arguments;
    try {
        arguments = std::make_unique<ProgramArguments>(argc, argv);
    } catch (InvalidArgumentsException &e) {
        die(e.what());
    }

    if (arguments->isHelpRequested())
        return handle_help_request(*arguments);

    std::string mode = arguments->getMode();
    if (programModes.find(mode) == programModes.end()) {
        std::cerr << "Unknown mode: " << mode << ". Type '" << arguments->getCmd() << " help' for help" << std::endl;
        return EXIT_FAILURE;
    }

    Parameters params = arguments->getParameters();
    ShapeFactory::initShapeClass(params.particleType, params.particleAttributes);

#ifdef _OPENMP
    omp_set_num_threads(params.ompThreads);
#endif

    auto &modeObject = *(programModes[mode]);
    modeObject.initializeForArguments(*arguments);
    modeObject.run();
    return EXIT_SUCCESS;
}
