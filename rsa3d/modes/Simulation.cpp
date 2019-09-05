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

        long page_size_kb = sysconf(_SC_PAGE_SIZE); // in case x86-64 is configured to use 2MB pages
        vm_usage = vsize / 1024.0 / 1024.0;
        resident_set = rss * page_size_kb / 1024.0 / 1024.0;
    }
}


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
    std::ofstream datFile(sFile);
    if (!datFile)
        die("Cannot open file " + sFile + " to store packing info");
    datFile.precision(std::numeric_limits<double>::digits10 + 1);

    std::size_t packingIndex{};
    Packing packing;
    do {
        std::size_t seed = params.from + packingIndex;
        packing = runSingleSimulation(seed, datFile);
    } while(this->continuePackingGeneration(packingIndex++, packing));

    this->postProcessSimulation();
}