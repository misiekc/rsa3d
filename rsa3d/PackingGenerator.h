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
#include "VoxelList.h"
#include "Surface.h"
#include <vector>
#include <map>

class PackingGenerator {
private:
	static double FACTOR_LIMIT;

	int seed;
	Parameters *params;
	std::vector<Shape *> packing;
	VoxelList *voxels;
	Surface *surface;

	int analyzeVoxels();
	int analyzeRegion(Voxel *v);
	bool isSaturated();
	double getFactor();
	void createPacking();

public:
	PackingGenerator(int seed, Parameters *params);
	virtual ~PackingGenerator();

	void run();
	std::vector<Shape *> * getPacking();
	void toPovray(std::string filename);
	static void toPovray(std::vector<Shape *> * packing, double size, std::string filename);
};

#endif /* PACKINGGENERATOR_H_ */
