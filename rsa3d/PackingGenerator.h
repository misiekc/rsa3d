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

    static void expandShapeOnBC(Packing *packing, const RSAShape *shape, double translation, size_t translateCoordIdx);

public:
    PackingGenerator(int seed, Parameters *params);

	virtual ~PackingGenerator();
    void run();

	Packing *getPacking();

	void testPacking(Packing *packing, double maxTime);
    void toFile(const std::string &filename);
	void restore(std::istream &f);

	/**
	 * @brieg Restores dynamically allocated packing from file and returns it
	 * @param filename file name to read from
	 * @return restored packing
	 */
	static Packing *fromFile(const std::string &filename);

    /**
     * @brief Stores a @a packing to given file. Dies if a file cannot be created.
     * @param packing a packing to store
     * @param filename the name of a file to store packing to
     */
    static void toFile(const Packing *packing, const std::string &filename);

	static void toPovray(Packing *packing, double size, VoxelList *voxels, const std::string &filename);
	static void toWolfram(Packing *packing, double size, VoxelList *voxels, const std::string &filename);

    /**
     * @brief Clones and translates shapes distant from the surface boundary not more than @a expandMargin x @a size
     * according to periodic boundary conditions
     * @param packing a packing to expand
     * @param size the size of the packing
     * @param expandMargin the distance (from 0 to 0.5) relative to @a size from the surface boundary from which shapes
     * will be translated
     */
    static void expandPackingOnPBC(Packing *packing, double size, double expandMargin);
};

#endif /* PACKINGGENERATOR_H_ */
