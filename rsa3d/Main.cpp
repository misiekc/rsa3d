#include "PackingGenerator.h"
#include "analizator/Analyzer.h"
#include "analizator/ExclusionZoneVisualizer.h"
#include "shape/ShapeFactory.h"
#include "utils/Quantity.h"
#include "utils/Assertions.h"
#include "ProgramArguments.h"


#include <unistd.h>
#include <sys/wait.h>
#include <cstring>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <memory>


namespace {
    //////////////////////////////////////////////////////////////////////////////
    //
    // process_mem_usage(double &, double &) - takes two doubles by reference,
    // attempts to read the system-dependent data for a process' virtual memory
    // size and resident set size, and return the results in KB.
    //
    // On failure, returns 0.0, 0.0
    void process_mem_usage(double &vm_usage, double &resident_set) {
        using std::ios_base;
        using std::ifstream;
        using std::string;

        vm_usage = 0.0;
        resident_set = 0.0;

        // 'file' stat seems to give the most reliable results
        //
        ifstream stat_stream("/proc/self/stat", ios_base::in);

        // dummy vars for leading entries in stat that we don't care about
        //
        string pid, comm, state, ppid, pgrp, session, tty_nr;
        string tpgid, flags, minflt, cminflt, majflt, cmajflt;
        string utime, stime, cutime, cstime, priority, nice;
        string O, itrealvalue, starttime;

        // the two fields we want
        //
        unsigned long vsize;
        long rss;

        stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                    >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                    >> utime >> stime >> cutime >> cstime >> priority >> nice
                    >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

        stat_stream.close();

        long page_size_kb = sysconf(_SC_PAGE_SIZE); // in case x86-64 is configured to use 2MB pages
        vm_usage = vsize / 1024.0 / 1024.0;
        resident_set = rss * page_size_kb / 1024.0 / 1024.0;
    }

    class Simulation {
    private:
        Packing runSingleSimulation(int seed, std::ofstream &dataFile);

    protected:
        Parameters params;

    public:
        explicit Simulation(Parameters params) : params{std::move(params)} { }

        void run();

        virtual bool continuePackingGeneration(std::size_t packingIndex, const Packing &packing) = 0;
        virtual void postProcessPacking(Packing &packing) = 0;
    };

    Packing Simulation::runSingleSimulation(int seed, std::ofstream &dataFile) {
        double vm, rss;

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        PackingGenerator pg(seed, &params);
        pg.run();
        Packing packing = pg.getPacking();
        this->postProcessPacking(packing);

        if (params.storePackings)
            packing.store(pg.getPackingFilename());

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        process_mem_usage(vm, rss);
        dataFile << seed << "\t" << packing.size() << "\t" << packing.back()->time << std::endl;
        dataFile.flush();

        auto generationSeconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
        std::cout << "[Simulation::runSingleSimulation] T: " << generationSeconds << "; VM: " << vm << "; RSS: " << rss;
        std::cout << std::endl << std::endl;

        return packing;
    }

    void Simulation::run() {
        std::string sFile = params.getPackingSignature() + ".dat";
        std::ofstream file(sFile);
        if (!file)
            die("Cannot open file " + sFile + " to store packing info");
        file.precision(std::numeric_limits<double>::digits10 + 1);

        std::size_t packingIndex{};
        Packing packing;
        do {
            std::size_t seed = params.from + packingIndex;
            packing = runSingleSimulation(seed, file);
        } while(this->continuePackingGeneration(packingIndex++, packing));
    }

    class DefaultSimulation : public Simulation {
    public:
        explicit DefaultSimulation(const Parameters &parameters) : Simulation(parameters) { }

        bool continuePackingGeneration(std::size_t packingIndex, const Packing &packing) override {
            return packingIndex < this->params.collectors - 1;
        }

        void postProcessPacking(Packing &packing) override { }
    };

    class BoundariesSimulation : public Simulation {
    private:
        std::size_t targerNumberOfParticles{};
        std::size_t particleCounter{};

    public:
        BoundariesSimulation(const Parameters &parameters, std::size_t targetNumberOfParticles)
                : Simulation(parameters), targerNumberOfParticles{targetNumberOfParticles}
        {
            Expects(targetNumberOfParticles > 0);
        }

        bool continuePackingGeneration(std::size_t packingIndex, const Packing &packing) override {
            return this->particleCounter < this->targerNumberOfParticles;
        }

        void postProcessPacking(Packing &packing) override {
            this->particleCounter += packing.size();
            std::cout << "[BoundariesSimulation::postProcessPacking] Fraction of particles generated: ";
            std::cout << static_cast<double>(this->particleCounter) / this->targerNumberOfParticles << std::endl;
        }
    };

    class AccuracySimulation : public Simulation {
    private:
        double targetAccuracy{};
        std::vector<double> packingFractions;
        Quantity meanPackingFraction;

    public:
        AccuracySimulation(const Parameters &parameters, double targetAccuracy) : Simulation(parameters),
                                                                                  targetAccuracy{targetAccuracy}
        {
            Expects(targetAccuracy > 0);
        }

        bool continuePackingGeneration(std::size_t packingIndex, const Packing &packing) override {
            return packingIndex < 2 || this->meanPackingFraction.error > this->targetAccuracy;
        }

        void postProcessPacking(Packing &packing) override {
            this->packingFractions.push_back(static_cast<double>(packing.size()) / this->params.sufraceVolume());
            this->meanPackingFraction.calculateFromSamples(this->packingFractions);
            std::cout << "[AccuracySimulation::postProcessPacking] Actual accuracy: ";
            std::cout << this->meanPackingFraction.error << std::endl;
            std::cout << "[AccuracySimulation::postProcessPacking] Target accuracy: ";
            std::cout << this->targetAccuracy << std::endl;
        }

        Quantity getMeanPackingFraction() const {
            return this->meanPackingFraction;
        }
    };

    class DensitySimulation : public Simulation {
    private:
        double targetPackingFraction{};

    public:
        DensitySimulation(const Parameters &parameters, double targetPackingFraction)
                : Simulation(parameters), targetPackingFraction{targetPackingFraction}
        {
            Expects(targetPackingFraction > 0.0 && targetPackingFraction <= 1.0);
        }

        bool continuePackingGeneration(std::size_t packingIndex, const Packing &packing) override {
            return packingIndex < this->params.collectors - 1;
        }

        void postProcessPacking(Packing &packing) override {
            double currentPackingFraction = packing.getParticlesVolume(this->params.surfaceDimension)
                                            / this->params.sufraceVolume();

            if (this->targetPackingFraction > currentPackingFraction) {
                std::cout << "[DensitySimulation::postProcessPacking] Target packing fraction: ";
                std::cout << this->targetPackingFraction << " > generated packing fraction: ";
                std::cout << currentPackingFraction << ". Keeping original packing." << std::endl;
            } else {
                std::cout << "[DensitySimulation::postProcessPacking] Initial packing fraction      : ";
                std::cout << currentPackingFraction << ", reducing... " << std::flush;
                double targetVolume = this->targetPackingFraction * this->params.sufraceVolume();
                packing.reducePackingVolume(targetVolume, this->params.surfaceDimension);
                std::cout << "done." << std::endl;

                double actualPackingFraction = packing.getParticlesVolume(this->params.surfaceDimension)
                                               / this->params.sufraceVolume();
                std::cout << "[DensitySimulation::postProcessPacking] Target packing fraction       : ";
                std::cout << this->targetPackingFraction << std::endl;
                std::cout << "[DensitySimulation::postProcessPacking] Actual packing fraction       : ";
                std::cout << actualPackingFraction << std::endl;
                std::cout << "[DensitySimulation::postProcessPacking] Last particle adsorption time : ";
                std::cout << packing.back()->time << std::endl;
            }
        }
    };
}

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

int simulate(const ProgramArguments &arguments) {
    if (!arguments.getPositionalArguments().empty())
        die(arguments.formatUsage(""));

    DefaultSimulation defaultSimulation(arguments.getParameters());
    defaultSimulation.run();

    /*std::string sFile = params.getPackingSignature() + ".dat";
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
    file.close();*/
    return 1;
}


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

void boundaries(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("<number of particles>"));

    std::size_t targetNumberOfParticles = std::stoul(positionalArguments[0]);
    Validate(targetNumberOfParticles > 0);

    BoundariesSimulation boundariesSimulation(arguments.getParameters(), targetNumberOfParticles);
    boundariesSimulation.run();
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

void accuracy(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 2)
        die(arguments.formatUsage("<accuracy> <out file>"));
    double accuracy = std::stod(positionalArguments[0]);
    Validate(accuracy > 0);

    Parameters params = arguments.getParameters();
    AccuracySimulation accuracySimulation(params, accuracy);
    accuracySimulation.run();
    Quantity meanPackingFraction = accuracySimulation.getMeanPackingFraction();

    std::string filename = positionalArguments[1];
    std::ofstream file(filename, std::ios_base::app);
    if (!file)
        die("Cannot open " + filename + " file to store packing fraction");
    file.precision(std::numeric_limits<double>::digits10 + 1);
    file << meanPackingFraction.value << "\t" << meanPackingFraction.error << "\t" << params.particleType;
    file << "\t" << params.particleAttributes << std::endl;
    file.close();
}

void density(const ProgramArguments &arguments) {
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        die(arguments.formatUsage("<target packing fraction>"));
    double targetPackingFraction = std::stod(positionalArguments[0]);
    Validate(targetPackingFraction > 0.0);
    Validate(targetPackingFraction <= 1.0);

    DensitySimulation densitySimulation(arguments.getParameters(), targetPackingFraction);
    densitySimulation.run();
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
    
    std::string mode = arguments->getMode();
    if (mode == "simulate")
        simulate(*arguments);
    else if (mode == "test")
        test(*arguments);
    else if (mode == "debug")
        debug(*arguments);
    else if (mode == "boundaries")
        boundaries(*arguments);
    else if (mode == "accuracy")
        accuracy(*arguments);
    else if (mode == "density")
        density(*arguments);
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

    return EXIT_SUCCESS;
}
