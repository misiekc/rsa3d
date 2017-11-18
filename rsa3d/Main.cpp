#include "PackingGenerator.h"
#include "Shape.h"
#include "shapes/Sphere.h"
#include "shapes/Cuboid.h"
#include "ShapeFactory.h"
#include "RND.h"
#include "analizator/Analyzer.h"
#include "tests/CuboidSpeedTest.h"
#include "tests/BoxFactory.h"
#include "tests/CuboidIntTest.h"
#include "tests/CuboidPointInsideTest.h"
#include "Config.h"

#include <unistd.h>
#include <sys/wait.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <memory>


void die(const std::string & reason)
{
    std::cerr << reason << std::endl;
    exit(EXIT_FAILURE);
}

void toFile(std::string filename, std::vector<Shape *> *packing) {
	std::ofstream file(filename, std::ios::binary);
	for (Shape *s : *packing) {
		s->store(file);
	}
	file.close();
}

std::vector<Shape *> * fromFile(unsigned char dim, std::string filename) {
	std::ifstream file(filename, std::ios::binary);
	std::vector<Shape *> * v = new std::vector<Shape *>;
	RND rnd(1);

	while (!file.eof()) {
		Shape *s = ShapeFactory::createShape(&rnd);
		s->restore(file);
		v->push_back(s);
	}
	v->pop_back();

	file.close();
	return v;
}

int simulate(Parameters *params) {
	PackingGenerator *pg;
	char buf[20];
	int pid = 0;
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, params->dimension));
	std::string size(buf);

	std::string sFile = "packing_" + params->particleType + "_"
			+ params->particleAttributes + "_" + size + ".dat";
	std::ofstream file(sFile);
	file.precision(std::numeric_limits<double>::digits10 + 1);
	std::vector<Shape *> *packing;

	for (int i = params->from; i < params->from + params->collectors; i++) {

		if (params->multiProcess){
			pid = fork();
			if (pid < 0){
				std::cout << "fork problem" << std::endl;
				i--;
				continue;
			}
			if (pid > 0)
				continue;
		}
		pg = new PackingGenerator(i, params);
		pg->run();
		packing = pg->getPacking();
		if (params->storePackings) {
			std::string sPackingFile = "packing_" + params->particleType + "_" + params->particleAttributes + "_" + size + "_" + std::to_string(i) + ".bin";
			toFile(sPackingFile, packing);
		}
		file << i << "\t" << packing->size() << "\t"
				<< (*packing)[packing->size() - 1]->time << std::endl;
		file.flush();
		delete pg;
		if (params->multiProcess && pid == 0)
			return 1;
	}
	if (params->multiProcess){
		for (int i = params->from; i < params->from + params->collectors; i++)
			wait(NULL);
	}
	file.close();
	return 1;
}

void boundaries(Parameters *params) {
	PackingGenerator *pg;
	char buf[20];
	unsigned long counter = 0UL, max = 200000000UL;
	int seed = 0;
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, params->dimension));
	std::string size(buf);
	std::string filename = "packing_" + params->particleType + "_" + params->particleAttributes + "_" + size + ".dat";
	std::ofstream file(filename);
	file.precision(std::numeric_limits<double>::digits10 + 1);
	std::vector<Shape *> *packing;
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


// Performs cuboid speedtests with parameters passed to process
//--------------------------------------------------------------------------------------------
int cube_speedtest_main(int argc, char **argv)
{    
    if (argc < 3)
        die("Usage: ./rsa cube_speedtest [input] [output = cube_speedtest.csv]");

    std::ifstream input(argv[2]);
    auto config = std::unique_ptr<Config> (Config::parse(input));
    input.close();    
    
    std::string output = "cube_speedtest.csv";
    if (argc == 4)
        output = argv[3];
    
    std::size_t pairs = config->getUnsignedInt("pairs");
    std::size_t repeats = config->getUnsignedInt("repeats");
    std::istringstream sizes(config->getString("box_sizes"));
    
	ShapeFactory::initShapeClass("Cuboid", "3 " + config->getString("cuboid_size"));
    BoxFactory * factory = BoxFactory::getInstance();
    std::vector<cube_speedtest::TestData> dataVector;
    double size;
    
    // Warm up and perform tests
    cube_speedtest::warmUp(factory);
    while (sizes >> size) {
        size /= 2;
        factory->setBoxSize(size, size, size);
        cube_speedtest::TestData data = cube_speedtest::perform(factory, pairs, repeats);
        dataVector.push_back(data);
    }
    
    // Print results
    std::cout << std::endl;
    std::cout << ">> Obtained results:" << std::endl << std::endl;
    for (auto data : dataVector) {
        data.printResults();
        std::cout << std::endl;
    }
    
    // Print aquired probabilities
    std::cout << ">> Probabilities for box sizes" << std::endl;
    std::cout << std::left;
    for (auto data : dataVector) {
        std::cout << std::setw(15) << data.overlapProb.value << " : " << data.factoryDesc << std::endl;
    }
    std::cout << std::endl;
    
    // Store to file
    std::cout << ">> Storing to file " << output << "..." << std::endl;
    std::ofstream file(output);
    cube_speedtest::TestData::toCsv(file, dataVector);
    file.close();
    
    return EXIT_SUCCESS;
}

// Performs cuboid intersection test with parameters passed to process
//--------------------------------------------------------------------------------------------
void cube_inttest_main(int argc, char ** argv)
{
    if (argc < 7)
        die("Usage: ./rsa cube_inttest [size_x] [size_y] [size_z] [box_halfsize] [max_tries]");
    
    double box_halfsize = std::stod(argv[5]);
    int max_tries = std::stoi(argv[6]);
    if (box_halfsize <= 0 || max_tries <= 0)
        die("Wrong input. Aborting.");
        
    std::stringstream init_stream;
    init_stream << "3 " << argv[2] << " " << argv[3] << " " << argv[4];
	ShapeFactory::initShapeClass("Cuboid", init_stream.str());
    BoxFactory * factory = BoxFactory::getInstance();
    factory->setBoxSize(box_halfsize, box_halfsize, box_halfsize);
    
    // Test SAT
    cube_inttest::Results results = cube_inttest::perform(factory, Cuboid::OverlapStrategy::SAT, max_tries);
    cube_inttest::print_results(results);
    std::cout << std::endl;
    results.free_missed_pairs();
    
    // Test TRI_TRI
    results = cube_inttest::perform(factory, Cuboid::OverlapStrategy::TRI_TRI, max_tries);
    cube_inttest::print_results(results);
    std::cout << std::endl;
    
    // Dump missed test to file
    if (results.missed > 0) {
        std::ofstream dump_file("inttest_dump.nb");
        if (dump_file) {
            cube_inttest::dump_missed_pairs(results, dump_file);
            results.free_missed_pairs();
            dump_file.close();
            std::cout << ">> Missed pairs dumped to inttest_dump.nb" << std::endl;
        } else {
            std::cout << ">> Could not write to inttest_dump.nb";
        }
    }
}

// Performs Cuboid::pointInside test
//--------------------------------------------------------------------------------------------
void cube_pitest_main(int argc, char ** argv)
{
    if (argc < 7)
        die("Usage: ./rsa cube_pi_test [size_x] [size_y] [size_z] [box_halfsize] [max_tries]");
    
    double box_halfsize = std::stod(argv[5]);
    int max_tries = std::stoi(argv[6]);
   
    if (box_halfsize <= 0 || max_tries <= 0)
        die("Wrong input. Aborting.");
        
	ShapeFactory::initShapeClass("Cuboid", std::string("3 ") + argv[2] + " " + argv[3] + " " + argv[4]);
    BoxFactory * factory = BoxFactory::getInstance();
    factory->setBoxSize(box_halfsize, box_halfsize, box_halfsize);
        
    cube_pitest::Results results = cube_pitest::perform(factory, max_tries);
    std::cout << std::endl;
    cube_pitest::print_results(results);
}

int main(int argc, char **argv) {
    if (argc < 2)
        die("No mode param. Aborting.");

    // Modes with custom environment
    if (strcmp(argv[1], "cube_speedtest") == 0) {
        cube_speedtest_main(argc, argv);
        return EXIT_SUCCESS;
    } else if (strcmp(argv[1], "cube_inttest") == 0) {
        cube_inttest_main(argc, argv);
        return EXIT_SUCCESS;
    } else if (strcmp(argv[1], "cube_pitest") == 0) {
        cube_pitest_main(argc, argv);
        return EXIT_SUCCESS;
    }
	
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
		std::vector<Shape *> *packing = fromFile(params.dimension, file);
		PackingGenerator::toPovray(packing, params.surfaceSize, NULL, file + ".pov");
		delete packing;
	}
	return 1;
}
