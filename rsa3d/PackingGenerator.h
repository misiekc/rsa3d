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

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
class PackingGenerator {
private:
	static double FACTOR_LIMIT;

	int seed;
	Parameters *params;
	std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> packing;
	VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *voxels;
	Surface *surface;

	int analyzeVoxels();
	int analyzeRegion(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v);
	void modifiedRSA(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s, Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v);
	bool isSaturated();
	double getFactor();
	void createPacking();

	void toPovray(const std::string &filename);
	void toWolfram(const std::string &filename);
	void toWolfram(double *da, const std::string &filename);

	void printRemainingVoxels(const std::string &prefix);


public:
	PackingGenerator(int seed, Parameters *params);
	virtual ~PackingGenerator();

	void run();
	std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> * getPacking();
	void toFile(const std::string &filename);
	static void toPovray(std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> * packing, double size, VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *voxels, const std::string &filename);
	static void toWolfram(std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> * packing, double size, VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *voxels, const std::string &filename);
};

#include "PackingGenerator.tpp"

#endif /* PACKINGGENERATOR_H_ */
