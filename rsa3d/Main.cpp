#include "PackingGenerator.h"
#include "analizator/Analyzer.h"
#include "ShapeFactory.h"

#include <unistd.h>
#include <sys/wait.h>
#include <cstring>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <memory>

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

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
   vm_usage     = vsize / 1024.0 / 1024.0;
   resident_set = rss * page_size_kb / 1024.0 / 1024.0;
}


void runSingleSimulation(int seed, Parameters *params, std::ofstream &dataFile){
	char buf[20];
	double vm, rss;
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, RSA_SPATIAL_DIMENSION));
	std::string size(buf);

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	PackingGenerator *pg = new PackingGenerator(seed, params);
	pg->run();
	Packing *packing = pg->getPacking();

	if (params->storePackings) {
		std::string sPackingFile = "packing_" + params->particleType + "_" + params->particleAttributes + "_" + size + "_" + std::to_string(seed) + ".bin";
		pg->toFile(sPackingFile);
	}
	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
	process_mem_usage(vm, rss);
	dataFile << seed << "\t" << packing->size() << "\t"	<< (*packing)[packing->size() - 1]->time << std::endl;
	dataFile.flush();
	std::cout << "T:" << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "; VM: " << vm << "; RSS: " << rss << std::endl;
	delete pg;
}

void debug(Parameters *params, char *cfile){
	std::string filename(cfile);
	char buf[20];
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, RSA_SPATIAL_DIMENSION));
	std::string size(buf);

	std::ifstream file(filename, std::ios::binary);
	if (!file)
		die("Cannot open file " + filename + " to restore packing generator");

	PackingGenerator *pg = new PackingGenerator(0, params);
	pg->restore(file);
	pg->run();
	delete pg;
}

int simulate(Parameters *params) {
	int pid = 0;
	char buf[20];
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, RSA_SPATIAL_DIMENSION));
	std::string size(buf);

	std::string sFile = "packing_" + params->particleType + "_" + params->particleAttributes + "_" + size + ".dat";
	std::ofstream file(sFile);
    if (!file)
        die("Cannot open file " + sFile + " to store packing info");
	file.precision(std::numeric_limits<double>::digits10 + 1);

	int seed = params->from;
	while(seed < params->from + params->collectors){
		int i;
		for(i=0; ( (i < params->generatorProcesses) && ((seed+i) < (params->from + params->collectors)) ); i++){
			pid = fork();
			if (pid < 0){
				std::cout << "fork problem" << std::endl;
				i--;
				continue;
			}
			if (pid==0){
				runSingleSimulation(seed + i, params, file);
				return 1;
			}
			if (pid > 0){
				continue;
			}
		}
		for(int j=0; j<i; j++){
			wait(NULL);
			seed++;
		}
	}
	file.close();
	return 1;
}


/**
 * @brief Creates packing until the total number of shapes (in all packings) does not exceed @a max
 */
void boundaries(Parameters *params, unsigned long max) {
	PackingGenerator *pg;
	char buf[20];
	unsigned long counter = 0;;
	int seed = 0;
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, RSA_SPATIAL_DIMENSION));
	std::string size(buf);
	std::string filename = "packing_" + params->particleType + "_" + params->particleAttributes + "_" + size + ".dat";
	std::ofstream file(filename);
	file.precision(std::numeric_limits<double>::digits10 + 1);
	Packing *packing;
	do {
		pg = new PackingGenerator(seed, params);
		pg->run();
		packing = pg->getPacking();
		file << seed << "\t" << packing->size() << "\t" << (*packing)[packing->size() - 1]->time << std::endl;
		file.flush();
		counter += pg->getPacking()->size();
		std::cout << (double) counter / (double) max << std::endl;
		seed++;
		delete pg;
	} while (counter < max);
	file.close();
}


int main(int argc, char **argv) {
    if (argc < 3)
        die("Usage: ./rsa <mode> <config> (additional parameters)");

    Parameters params(argv[2]);
    ShapeFactory::initShapeClass(params.particleType, params.particleAttributes);

    std::string mode(argv[1]);
    if (mode == "simulate") {
        simulate(&params);
    } else if (mode == "test") {
        if (argc < 5)   die("Usage: ./rsa test <config> <file in> <max time>");
        std::string filename(argv[3]);
        Packing *packing = PackingGenerator::fromFile(filename);
        PackingGenerator *pg = new PackingGenerator(1, &params);
        pg->testPacking(packing, atof(argv[4]));
        delete pg;
        delete packing;
    } else if (mode == "debug") {
        debug(&params, argv[3]);
    } else if (mode == "boundaries") {
        boundaries(&params, atof(argv[3]));
    } else if (mode == "analyze") {
        Analyzer an(&params);
        an.analyzePackingsInDirectory(argv[3], 0.01, 1.0);
    } else if (mode == "povray") {
        std::string file(argv[3]);
        Packing *packing = PackingGenerator::fromFile(argv[3]);
        PackingGenerator::toPovray(packing, params.surfaceSize, nullptr, file + ".pov");
        delete packing;
    } else if (mode == "wolfram") {
        std::string file(argv[3]);
        Packing *packing = PackingGenerator::fromFile(file);
        PackingGenerator::toWolfram(packing, params.surfaceSize, nullptr, file + ".nb");
        delete packing;
    } else if (mode == "bc_expand") {
        if (argc < 4)   die("Usage: ./rsa bc_expand <config> <file in> (file out = file in)");
        std::string fileIn(argv[3]);
        std::string fileOut = (argc < 5 ? fileIn : argv[4]);
        auto packing = PackingGenerator::fromFile(fileIn);
        PackingGenerator::expandPackingOnPBC(packing, params.surfaceSize, 0.1);
        PackingGenerator::toFile(packing, fileOut);
        delete packing;
    } else if (mode == "overlap") {
        if (argc < 4)   die("Usage: ./rsa overlap <config> <file in>");
        std::string fileIn(argv[3]);
        auto packing = PackingGenerator::fromFile(fileIn);
        PackingGenerator::testPackingOverlaps(packing);
        delete packing;
    } else {
		std::cerr << "Unknown mode: " << argv[1] << std::endl;
        return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
