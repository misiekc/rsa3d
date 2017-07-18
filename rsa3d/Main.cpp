#include "PackingGenerator.h"
#include "Shape.h"
#include "shapes/Sphere.h"
#include "shapes/Cuboid.h"
#include "ShapeFactory.h"
#include "RND.h"
#include "analizator/Analyzer.h"

#include <string.h>
#include <fstream>
#include <iostream>

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
	std::sprintf(buf, "%.0f", pow(params->surfaceSize, params->dimension));
	std::string size(buf);

	std::string sFile = "packing_" + params->particleType + "_"
			+ params->particleAttributes + "_" + size + ".dat";
	std::ofstream file(sFile);
	file.precision(std::numeric_limits<double>::digits10 + 1);
	std::vector<Shape *> *packing;

	for (int i = params->from; i < params->from + params->collectors; i++) {
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

int main(int argc, char **argv) {
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
