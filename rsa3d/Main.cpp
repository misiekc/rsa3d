#include "PackingGenerator.h"
#include "analizator/Analyzer.h"
#include "analizator/ExclusionZoneVisualizer.h"
#include "shape/ShapeFactory.h"
#include "Utils.h"


#include <unistd.h>
#include <sys/wait.h>
#include <cstring>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <memory>
#include <dirent.h>

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

void makeDatFileForPackingsInDirectory(Parameters *params, char *sdir){
	char prefix[] = "packing";

	char buf[20];
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, params->surfaceDimension));
	std::string size(buf);
	std::string sFile;

	if (params->particleAttributes.length()<100)
	    sFile = "packing_" + params->particleType + "_" + replaceAll(params->particleAttributes, " ", "_") + "_" + size + ".dat";
	else
        sFile = "packing_" + params->particleType + "_" + replaceAll(params->particleAttributes, " ", "_").substr(0, 100) + "_" + size + ".dat";
	std::ofstream dataFile(sFile);
	if (!dataFile)
    	die("Cannot open file " + sFile + " to store packing info");
	dataFile.precision(std::numeric_limits<double>::digits10 + 1);

	DIR *dir = opendir(sdir);
	dirent *de;
	std::string dirname(sdir);
	while ((de = readdir(dir)) != nullptr){
		if (std::strncmp(prefix, de->d_name, strlen(prefix))==0){
			int no1 = lastIndexOf(de->d_name, '_');
			int no2 = lastIndexOf(de->d_name, '.');
			char seed[10];
			strncpy(seed, de->d_name+no1+1, no2-no1-1);
			seed[no2-no1-1] = '\0';
			std::string filename(de->d_name);
			Packing packing;
			packing.restore(dirname + "/" + filename);
			dataFile << seed << "\t" << packing.back()->no << "\t"	<< packing.back()->time  <<  std::endl;
			std::cout << ".";
			std::cout.flush();
		}
	}
	(void)closedir(dir);
	std::cout << std::endl;
}

void runSingleSimulation(int seed, Parameters *params, std::ofstream &dataFile){
	char buf[20];
	double vm, rss;
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, params->surfaceDimension));
	std::string size(buf);

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	PackingGenerator pg(seed, params);
	pg.run();
	const Packing &packing = pg.getPacking();

	if (params->storePackings) {
        std::string sPackingFile;
        if (params->particleAttributes.length()<100)
    		sPackingFile = "packing_" + params->particleType + "_" + replaceAll(params->particleAttributes, " ", "_") + "_" + size + "_" + std::to_string(seed) + ".bin";
        else
            sPackingFile = "packing_" + params->particleType + "_" + replaceAll(params->particleAttributes, " ", "_").substr(0, 100) + "_" + size + "_" + std::to_string(seed) + ".bin";
		pg.getPacking().store(sPackingFile);
	}
	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
	process_mem_usage(vm, rss);
	dataFile << seed << "\t" << packing.size() << "\t"	<< packing.back()->time << std::endl;
	dataFile.flush();
	std::cout << "T:" << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "; VM: " << vm << "; RSS: " << rss << std::endl;
}

void debug(Parameters *params, char *cfile){
	std::string filename(cfile);
	char buf[20];
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, params->surfaceDimension));
	std::string size(buf);

	std::ifstream file(filename, std::ios::binary);
	if (!file)
		die("Cannot open file " + filename + " to restore packing generator");

	PackingGenerator pg(0, params);
	pg.restore(file);
	pg.run();
}

int simulate(Parameters *params) {
	int pid = 0;
	char buf[20];
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, params->surfaceDimension));
	std::string size(buf);


    std::string sFile;
    if (params->particleAttributes.length()<100)
        sFile = "packing_" + params->particleType + "_" + replaceAll(params->particleAttributes, " ", "_") + "_" + size + ".dat";
    else
        sFile = "packing_" + params->particleType + "_" + replaceAll(params->particleAttributes, " ", "_").substr(0, 100) + "_" + size + ".dat";
    std::ofstream file(sFile);
	if (!file)
    	die("Cannot open file " + sFile + " to store packing info");
	file.precision(std::numeric_limits<double>::digits10 + 1);

	std::size_t seed = params->from;

	if (params->generatorProcesses>1){
		while(seed < params->from + params->collectors){
			std::size_t i;
			for(i=0; ( (i < params->generatorProcesses) && ((seed+i) < (params->from + params->collectors)) ); i++){
				pid = fork();
				if (pid < 0){
					std::cout << "fork problem" << std::endl;
					i--;
					continue;
				}
				if (pid==0){
					runSingleSimulation(static_cast<int>(seed + i), params, file);
					return 1;
				}
				if (pid > 0){
					continue;
				}
			}
			for(std::size_t j=0; j<i; j++){
				wait(nullptr);
				seed++;
			}
		}
	}else{
		for(std::size_t i=0; i<params->collectors; i++){
			runSingleSimulation(static_cast<int>(params->from + i), params, file);
		}
	}
	file.close();
	return 1;
}


/**
 * @brief Creates packing until the total number of shapes (in all packings) does not exceed @a max
 */
void boundaries(Parameters *params, unsigned long max) {
	char buf[20];
    unsigned long counter = 0;;
    int seed = 0;
    std::sprintf(buf, "%.0f", pow(params->surfaceSize, params->surfaceDimension));
    std::string size(buf);
    std::string sFile;
    if (params->particleAttributes.length()<100)
        sFile = "packing_" + params->particleType + "_" + replaceAll(params->particleAttributes, " ", "_") + "_" + size + ".dat";
    else
        sFile = "packing_" + params->particleType + "_" + replaceAll(params->particleAttributes, " ", "_").substr(0, 100) + "_" + size + ".dat";
    std::ofstream file(sFile);
    file.precision(std::numeric_limits<double>::digits10 + 1);
    do {
        PackingGenerator pg(seed, params);
        pg.run();
        const Packing &packing = pg.getPacking();
		file << seed << "\t" << packing.size() << "\t" << packing.back()->time << std::endl;
		file.flush();
		counter += packing.size();
		std::cout << (double) counter / (double) max << std::endl;
		seed++;
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
        if (argc < 5)   die("Usage: ./rsa test <input> <file in> <max time>");
        std::string filename(argv[3]);
        Packing packing;
        packing.restore(filename);
        PackingGenerator pg(1, &params);
        pg.testPacking(packing, std::stod(argv[4]));
    } else if (mode == "debug") {
        debug(&params, argv[3]);
    } else if (mode == "boundaries") {
        boundaries(&params, std::stoul(argv[3]));
    } else if (mode == "dat") {
        makeDatFileForPackingsInDirectory(&params, argv[3]);
    } else if (mode == "analyze") {
    	if (argc < 4) die("Usage: ./rsa analyze <input> <directory> (correlations range = 10; 0 - no corr output)");
    	double corrRange = (argc >= 5) ? std::stod(argv[4]) : 10.0;
    	Validate(corrRange >= 0);
        Analyzer an(&params);
		an.analyzePackingsInDirectory(argv[3], 0.01, 1.0, corrRange);
    } else if (mode == "povray") {
        if (argc < 4)   die("Usage: ./rsa povray <input> <file in>");
        std::string file(argv[3]);
        Packing packing;
        packing.restore(file);
        PackingGenerator::toPovray(packing, params.surfaceSize, nullptr, file + ".pov");
    } else if (mode == "wolfram") {
        if (argc < 4)   die("Usage: ./rsa wolfram <input> <file in>");
        std::string file(argv[3]);
        Packing packing;
        packing.restore(file);
        PackingGenerator::toWolfram(packing, params.surfaceSize, nullptr, file + ".nb");
    } else if (mode == "bc_expand") {
        if (argc < 4)   die("Usage: ./rsa bc_expand <config> <file in> (file out = file in)");
        std::string fileIn(argv[3]);
        std::string fileOut = (argc < 5 ? fileIn : argv[4]);
        Packing packing;
        packing.restore(fileIn);
        packing.expandOnPBC(params.surfaceSize, 0.1);
        packing.store(fileOut);
    } else if (mode == "exclusion_zones"){
        if (argc < 5)   die("Usage: ./rsa test <input> <packign file> <output file>");
    	std::string packingFile(argv[3]);
    	std::string outputFile(argv[4]);
    	ExclusionZoneVisualizer::main(params, packingFile, outputFile);
    } else {
		std::cerr << "Unknown mode: " << argv[1] << std::endl;
        return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
