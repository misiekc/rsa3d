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

template <ushort DIMENSION>
class PackingGenerator {
private:
	static double FACTOR_LIMIT;

	int seed;
	Parameters *params;
	std::vector<Shape<DIMENSION> *> packing;
	VoxelList<DIMENSION> *voxels;
	Surface *surface;

	int analyzeVoxels();
	int analyzeRegion(Voxel<DIMENSION, ANGULAR_DIMENSION> *v);
	void modifiedRSA(Shape<DIMENSION> *s, Voxel<DIMENSION, ANGULAR_DIMENSION> *v);
	bool isSaturated();
	double getFactor();
	void createPacking();

public:
	PackingGenerator(int seed, Parameters *params);
	virtual ~PackingGenerator();

	void run();
	std::vector<Shape<DIMENSION> *> * getPacking();
	void toPovray(std::string filename);
	static void toPovray(std::vector<Shape<DIMENSION> *> * packing, double size, VoxelList<DIMENSION> *voxels, std::string filename);
};

#include "PackingGenerator.tpp"

#endif /* PACKINGGENERATOR_H_ */
