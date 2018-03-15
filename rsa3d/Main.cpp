#include "PackingGenerator.h"
#include "analizator/Analyzer.h"
#include "tests/CuboidSpeedTest.h"
#include "tests/CuboidPointInsideTest.h"
#include "tests/CuboidIntTest.h"
#include "tests/AnisotropicShape2DExclusionTest.h"
#include "tests/VectorSpeedTest.h"

#include <unistd.h>
#include <sys/wait.h>
#include <cstring>
#include <iomanip>


void toFile(const std::string &filename, std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *packing) {
	std::ofstream file(filename, std::ios::binary);
	for (Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s : *packing) {
		s->store(file);
	}
	file.close();
}

std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> * fromFile(const std::string &filename) {
	std::ifstream file(filename, std::ios::binary);
	if (!file)
		die("Cannot open file " + filename + " to restore packing");

	auto * v = new std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *>;
	RND rnd(1);

	while (!file.eof()) {
		Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s = ShapeFactory::createShape(&rnd);
		s->restore(file);
		v->push_back(s);
	}
	v->pop_back();

	file.close();
	return v;
}

void runSingleSimulation(int seed, Parameters *params, std::ofstream &dataFile){
	char buf[20];
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, RSA_SPATIAL_DIMENSION));
	std::string size(buf);

	auto *pg = new PackingGenerator<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>(seed, params);
	pg->run();
	std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *packing = pg->getPacking();

	if (params->storePackings) {
		std::string sPackingFile = "packing_" + params->particleType + "_" + params->particleAttributes + "_" + size + "_" + std::to_string(seed) + ".bin";
		toFile(sPackingFile, packing);
	}
	dataFile << seed << "\t" << packing->size() << "\t"	<< (*packing)[packing->size() - 1]->time << std::endl;
	dataFile.flush();
	delete pg;
}

int simulate(Parameters *params) {
	int pid = 0;
	char buf[20];
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, RSA_SPATIAL_DIMENSION));
	std::string size(buf);

	std::string sFile = "packing_" + params->particleType + "_" + params->particleAttributes + "_" + size + ".dat";
	std::ofstream file(sFile);
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

void boundaries(Parameters *params) {
	PackingGenerator<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *pg;
	char buf[20];
	unsigned long counter = 0UL, max = 200000000UL;
	int seed = 0;
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, RSA_SPATIAL_DIMENSION));
	std::string size(buf);
	std::string filename = "packing_" + params->particleType + "_" + params->particleAttributes + "_" + size + ".dat";
	std::ofstream file(filename);
	file.precision(std::numeric_limits<double>::digits10 + 1);
	std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *packing;
	do {
		pg = new PackingGenerator<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>(seed, params);
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
    if (argc < 2)
        die("No mode param. Aborting.");

    // Modes with custom environment
    if (strcmp(argv[1], "cube_speedtest") == 0)
		return cube_speedtest::main(argc, argv);
    else if (strcmp(argv[1], "cube_inttest") == 0)
		return cube_inttest::main(argc, argv);
    else if (strcmp(argv[1], "cube_pitest") == 0)
		return cube_pitest::main(argc, argv);
    else if (strcmp(argv[1], "as2d_extest") == 0)
        return as2d_extest::main(argc, argv);
    else if (strcmp(argv[1], "vec_speedtest") == 0)
        return vec_speedtest::main(argc, argv);
	
    // Modes with standard environment
    if (argc < 3)
        die("No input file given. Aborting.");
    
	Parameters params(argv[2]);
	ShapeFactory::initShapeClass(params.particleType, params.particleAttributes);

	if (strcmp(argv[1], "simulate")==0) {
		simulate(&params);
	} else if (strcmp(argv[1], "boundaries")==0) {
		boundaries(&params);
	} else if (strcmp(argv[1], "analyze")==0) {
		Analyzer an(&params);
		an.analyzePackingsInDirectory(argv[3], 0.01, 1.0);
	} else if (strcmp(argv[1], "povray")==0) {
		std::string file(argv[3]);
		std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *packing = fromFile(file);
		PackingGenerator<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>::toPovray(packing, params.surfaceSize, nullptr, file + ".pov");
		delete packing;
	} else if (strcmp(argv[1], "wolfram")==0) {
		std::string file(argv[3]);
		std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *packing = fromFile(file);
		PackingGenerator<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>::toWolfram(packing, params.surfaceSize, nullptr, file + ".nb");
		delete packing;
	}else{
		std::cerr << "Unknown mode: " << argv[1] << std::endl;
	}
	return 1;
}
