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
	std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> packing;
	VoxelList *voxels;
	Surface *surface;

	double spatialSize;
	double angularSize;

	int analyzeVoxels(unsigned short depth);
	int analyzeRegion(Voxel *v);
	void modifiedRSA(Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s, Voxel *v);
	bool isSaturated();
	double getFactor();
	bool isInside(double *position, double *orientation);
	void createPacking();

	void toPovray(const std::string &filename);
	void toWolfram(const std::string &filename);
	void toWolfram(double *da, const std::string &filename);

	void printRemainingVoxels(const std::string &prefix);

	void store(std::ostream &f) const;

public:
	PackingGenerator(int seed, Parameters *params);
	virtual ~PackingGenerator();

	void run();
	std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> * getPacking();

	void testPacking(std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *packing, double maxTime);

	void toFile(const std::string &filename);
	static void toPovray(std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> * packing, double size, VoxelList *voxels, const std::string &filename);
	static void toWolfram(std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> * packing, double size, VoxelList *voxels, const std::string &filename);

	void restore(std::istream &f);
};

#endif /* PACKINGGENERATOR_H_ */
