/*
 * PackingGenerator.h
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#ifndef PACKINGGENERATOR_H_
#define PACKINGGENERATOR_H_

#include "RND.h"
#include "Shape.h"
#include "BoundaryConditions.h"
#include "Parameters.h"
#include <vector>
#include <map>

class PackingGenerator {
private:
	int seed;
	Parameters *params;
	std::vector<Shape *> packing;
	void createPacking();

public:
	PackingGenerator(int seed, Parameters *params);
	virtual ~PackingGenerator();
	void run();
	std::vector<Shape *> * getPacking();
};

#endif /* PACKINGGENERATOR_H_ */
