/*
 * VoxelList.h
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#ifndef VOXELLIST_H_
#define VOXELLIST_H_

#include "Voxel.h"
#include "shape/Shape.h"
#include "RND.h"
#include "BoundaryConditions.h"
#include "Positioned.h"
#include <vector>

#include "NeighbourGrid.h"

class VoxelList {

private:
	// to optimize memory occupations voxels are initialized in the first invoke split()
	bool voxelsInitialized;

	// allows voxels to overlap - only for testing purposes and normally should be set to 1
	const double dxFactor = 1.0; // 1.0000000001;

	// neighbour grid structure for voxels. Needed to quickly find a voxel using its location
	NeighbourGrid<Voxel>* voxelNeighbourGrid;

	// probability distrubutions for drawing position inside a voxel
	std::uniform_real_distribution<double> *spatialDistribution;
	std::uniform_real_distribution<double> *angularDistribution;

	// array of top level voxel. If a top level voxel becomes inactive (due to shape placement, all its child voxels becomes obsolete
	bool* activeTopLevelVoxels;

	// returns initial size of a vovel. It is not grater than d and be an integer power of 2 (due to numerical issues)
	double findFloorSize(double d);

	// returns initial size of a vovel. It is not smaller than d and be an integer power of 2 (due to numerical issues)
	double findCeilSize(double d);

	// returns number of elements of size cellSize needed to cover a desired range
	size_t findArraySize(double range, double cellSize);

	// fills neigbour grid with voxels
	void rebuildNeighbourGrid();

	// for a given voxel returns index of its root
	int getIndexOfTopLevelVoxel(const RSAVector &da);

	// initialize voxels, returns number of initial voxels
	unsigned int initVoxels(RSABoundaryConditions *bc, NeighbourGrid<const RSAShape> *nl);

	// checks consistency of indexes of root voxels
	void checkTopLevelVoxels();


	// voxels array will not have NULLs between pointers to objects - they can appear when splitting or analyze voxels in parallel
	// returns number of not null values in list
	size_t compactVoxelArray(Voxel **list, size_t endIndex);

	// finds voxel containing given point - not used - only for debugging
	Voxel * findVoxel(Voxel **list, size_t listSize, const RSAVector &pos,
					  const RSAOrientation &angle);

protected:
	// disables voxel list for debug purposes - typically is false
	bool disabled;

	int surfaceDimension;
	Voxel** voxels;
	size_t length;

	double initialVoxelSize;
	double initialAngularVoxelSize;
	double spatialVoxelSize;
	double angularVoxelSize;


	double angularRange;
	double spatialRange;

	size_t beginningVoxelNumber;


	// checks if a top level voxel for voxel v is active (if not v should be removed
	bool isTopLevelVoxelActive(Voxel *v);

	bool isVoxelInsidePacking(Voxel *v);
	virtual bool isVoxelInsideExclusionZone(Voxel *v, double spatialSize, double angularSize,
                                    std::vector<const RSAShape*> *shapes, RSABoundaryConditions *bc,
                                    unsigned short depth = 0);
	bool isVoxelInsideExclusionZoneOld(Voxel *v, double spatialSize, double angularSize,
                                    std::vector<const RSAShape*> *shapes, RSABoundaryConditions *bc,
                                    unsigned short depth = 0);

	void splitVoxel(Voxel *v, double spatialSize, double angularSize, Voxel **vRes);

	virtual bool analyzeVoxel(Voxel *v, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc, double spatialSize, double angularSize, unsigned short depth = 0);

public:

	/**
	 * @brief Constructor
	 * @param packingSpatialSize packing size to be covered by voxels
	 * @param requestedSpatialVoxelSize suggested initial size of a voxel. Initial size of a allocated voxels will not be larger than the requested one.
	 * @param shapeAngularRange typpically 2*M_PI. Can be smaller for shapes with rotational symmetry. For example for squares it should be M_PI/2.0, and for ellipses, sherocylinders or rectangles M_PI.
	 * @param requestedAngularVoxelSize suggested initial angular size of a voxel. Initial angular size of allocaced voxels will not be larger than requested one.
	 */
	VoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange, double requestedAngularVoxelSize);

	void disable();

	virtual ~VoxelList();

	void getNeighbours(std::vector<Voxel *> *result, const RSAVector &da);
	void removeTopLevelVoxel(Voxel *v);

	size_t analyzeVoxels(RSABoundaryConditions *bc, NeighbourGrid<const RSAShape> *nl, unsigned short depth);


	double splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc);

	Voxel *getRandomVoxel(RND *rnd);
	Voxel *getVoxel(int i);
	Voxel *getVoxel(const RSAVector &pos, const RSAOrientation &angle);
	void getRandomPositionAndOrientation(RSAVector *position, RSAOrientation *orientation, Voxel *v, RND *rnd);
	double getVoxelSize();
	double getVoxelAngularSize();
	Voxel* get(int i);
	size_t getLength() const;
	double getVoxelsSurface();
	std::string toPovray();
	std::string toWolfram();

	void store(std::ostream &f) const;
	void restore(std::istream &f);

};


#endif /* VOXELLIST_H_ */
