//
// Created by pkua on 05.09.2019.
//

#include <chrono>
#include <fstream>
#include <unistd.h>

#include "Simulation.h"
#include "../PackingGenerator.h"


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
        #ifdef _SC_PAGE_SIZE
            long page_size_kb = sysconf(_SC_PAGE_SIZE); // in case x86-64 is configured to use 2MB pages
        #else
            long page_size_kb = 2048;
        #endif

        vm_usage = vsize / 1024.0 / 1024.0;
        resident_set = rss * page_size_kb / 1024.0 / 1024.0;
    }
}


Packing Simulation::runSingleSimulation(unsigned int seed, std::size_t collector, std::ofstream &dataFile) {
    double vm, rss;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    PackingGenerator pg(seed, collector, &params);
    Packing packing;
    if (access(pg.getPackingFilename().c_str(), F_OK)!=0){
    	pg.run();
    	packing = pg.getPacking();
    	this->postProcessPacking(packing);

    	if (this->params.storePackings)
    		packing.store(pg.getPackingFilename());
        dataFile << collector << "\t" << packing.size() << "\t" << packing.back()->time << std::endl;
        dataFile.flush();
    }else{
    	packing.restore(pg.getPackingFilename());
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    process_mem_usage(vm, rss);

    auto generationSeconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
    std::cout << "[Simulation::runSingleSimulation] T: " << generationSeconds << "; VM: " << vm << "; RSS: " << rss;
    std::cout << std::endl << std::endl;

    return packing;
}

void Simulation::run() {
    std::string datFilename = params.getPackingSignature() + ".dat";
    std::ofstream datFile;
    if (access(datFilename.c_str(), F_OK)==0)
        datFile.open(datFilename, std::ios_base::app);
    else
        datFile.open(datFilename, std::ios_base::out);
    if (!datFile)
        die("Cannot open file " + datFilename + " to store packing info");
    datFile.precision(std::numeric_limits<double>::digits10 + 1);

    std::size_t packingIndex{};
    Packing packing;
    unsigned long seedOrigin = this->getSeedOrigin();
    std::cout << "[Simulation::run] Using seed origin: " << seedOrigin << std::endl;
    do {
        std::size_t collector = params.from + packingIndex;
        unsigned long seed = seedOrigin + collector;
        packing = runSingleSimulation(seed, collector, datFile);
    } while(this->continuePackingGeneration(packingIndex++, packing));

    this->postProcessSimulation();
}

unsigned long Simulation::getSeedOrigin() const {
    if (this->params.seedOrigin == "random")
        return std::random_device{}();
    else
        return std::stoul(this->params.seedOrigin);
}

void Simulation::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << "This is one of simulation modes. Each one generates RSA packings based on given parameters." << std::endl;
    out << "It outputs a file" << std::endl;
    out << "    packing_[particle name]_[particle attributes]_[surface volume].dat" << std::endl;
    out << "whose rows correspond to subsequent packings and columns are: number of collector, number of" << std::endl;
    out << "particles and time of adding the last particle. Moreover, if specified, the packings are" << std::endl;
    out << "stored as" << std::endl;
    out << "    packing_[particle name]_[particle attributes]_[surface volume]_[number of collector].bin" << std::endl;
    out << std::endl;
    out << "The stopping conditions and additional operations specific to this mode:" << std::endl;
}
