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

/**
 * @brief A short form of vector of pointers to Shape with current parameters representing a packing
 */
using Packing = std::vector<RSAShape*>;

class PackingGenerator {
private:
	static double FACTOR_LIMIT;

	int seed;
	Parameters *params;
	Packing packing;
	VoxelList *voxels;
	Surface *surface;

	double spatialSize;
	double angularSize;

	int analyzeVoxels(unsigned short depth);
	int analyzeRegion(Voxel *v);
	void modifiedRSA(RSAShape *s, Voxel *v);
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
	Packing *getPacking();

	void testPacking(Packing *packing, double maxTime);

	void toFile(const std::string &filename);
	static void toPovray(Packing *packing, double size, VoxelList *voxels, const std::string &filename);
	static void toWolfram(Packing *packing, double size, VoxelList *voxels, const std::string &filename);

	void restore(std::istream &f);
};

#endif /* PACKINGGENERATOR_H_ */
