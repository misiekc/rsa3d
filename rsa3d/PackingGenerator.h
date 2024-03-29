/*
 * PackingGenerator.h
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#ifndef PACKINGGENERATOR_H_
#define PACKINGGENERATOR_H_

#include "RND.h"
#include "shape/Shape.h"
#include "BoundaryConditions.h"
#include "Parameters.h"
#include "VoxelList.h"
#include "Surface.h"
#include "Packing.h"
#include <vector>
#include <map>


class PackingGenerator {
private:
	static double FACTOR_LIMIT;

	int seed;
	std::size_t collector{};
	Parameters params;
	Packing packing;
	VoxelList *voxels;
	Surface *surface;

	double spatialSize;
	double angularSize;

    void modifiedRSA(RSAShape *s, Voxel *v);
	bool isSaturated();
	double getFactor();
	bool isInside(const RSAVector &position, RSAOrientation &orientation);
	void createPacking(Packing *packing);

	void toPovray(const std::string &filename);
    void toWolfram(const std::string &filename);
	void toWolfram(const RSAVector &da, const std::string &filename);

	void printRemainingVoxels(const std::string &prefix);

	void store(std::ostream &f) const;

public:
    PackingGenerator(int seed, std::size_t collector, const Parameters *params);

	virtual ~PackingGenerator();
    void run(Packing *packing=nullptr);

	const Packing &getPacking();

	void testPacking(const Packing &packing, double maxTime);
	void restore(std::istream &f);
	std::string getPackingFilename() const;

	static void toPovray(Packing packing, double size, VoxelList *voxels, bool drawPBC,
                         const std::string &filename);
    static void toWolfram(Packing packing, double size, VoxelList *voxels, bool isPeriodicImage,
                          const std::string &filename);
	static std::vector<std::string> findPackingsInDir(const std::string &dirName);

    bool generationCompleted(size_t missCounter, double t);

    [[nodiscard]] const Surface &getSurface() const { return *this->surface; }
};

#endif /* PACKINGGENERATOR_H_ */
