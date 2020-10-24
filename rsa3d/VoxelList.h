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

	// allows voxels to overlap - only for testing purposes and normally should be set to 1
	const double dxFactor = 1.0; // 1.0000000001;

	// array of top level voxel. If a top level voxel becomes inactive (due to shape placement, all its child voxels becomes obsolete
	bool* activeTopLevelVoxels;

	double spatialVoxelSize;
	double angularVoxelSize;

	// probability distrubutions for drawing position inside a voxel
	std::uniform_real_distribution<double> *spatialDistribution;
	std::uniform_real_distribution<double> *angularDistribution;

	// returns initial size of a vovel. It is not grater than d and be an integer power of 2 (due to numerical issues)
	double findFloorSize(double d);

	// returns initial size of a vovel. It is not smaller than d and be an integer power of 2 (due to numerical issues)
	double findCeilSize(double d);

	// returns number of elements of size cellSize needed to cover a desired range
	size_t findArraySize(double range, double cellSize);

	// for a given voxel returns index of its root
	int getIndexOfTopLevelVoxel(const RSAVector &da);

	// checks consistency of indexes of root voxels
	void checkTopLevelVoxels();

	// sets status basing on existing voxels
	void refreshTopLevelVoxels();

	// finds voxel containing given point - not used - only for debugging
	Voxel * findVoxel(Voxel **list, size_t listSize, const RSAVector &pos,
					  const RSAOrientation &angle);
protected:
	// disables voxel list for debug purposes - typically is false
	bool disabled;
	// to optimize memory occupations voxels are initialized in the first invoke split()
	bool voxelsInitialized;

	int surfaceDimension;
	Voxel** voxels;
	size_t length;

	double initialVoxelSize;
	double initialAngularVoxelSize;

	double angularRange;
	double spatialRange;

	size_t beginningVoxelNumber;

	// neighbour grid structure for voxels. Needed to quickly find a voxel using its location
	NeighbourGrid<Voxel>* voxelNeighbourGrid;

	virtual void allocateVoxels(size_t size);

	// initialize voxels, returns number of initial voxels
	unsigned int initVoxelsOld(RSABoundaryConditions *bc, NeighbourGrid<const RSAShape> *nl);
	unsigned int initVoxels(RSABoundaryConditions *bc, NeighbourGrid<const RSAShape> *nl);

	// checks if a top level voxel for voxel v is active (if not v should be removed
	bool isTopLevelVoxelActive(Voxel *v);

	virtual bool isVoxelInsidePacking(const Voxel *v, double spatialSize) const;
	virtual bool isVoxelInsideExclusionZone(Voxel *v, double spatialSize, double angularSize,
                                    std::vector<const RSAShape*> *shapes, RSABoundaryConditions *bc,
                                    unsigned short depth = 0);
	bool isVoxelInsideExclusionZoneOld(Voxel *v, double spatialSize, double angularSize,
                                    std::vector<const RSAShape*> *shapes, RSABoundaryConditions *bc,
                                    unsigned short depth = 0);

	void splitVoxel(Voxel *v, double spatialSize, double angularSize, Voxel **vRes);

	// fills neigbour grid with voxels
	void rebuildNeighbourGrid();

	bool analyzeVoxel(Voxel *v, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc, double spatialSize, double angularSize, unsigned short depth = 0);

	virtual void moveVoxelInList(size_t from, size_t to);

	// voxels array will not have NULLs between pointers to objects - they can appear when splitting or analyze voxels in parallel
	// returns number of not null values in list
	void compactVoxelArray();

public:

	/**
	 * constants returned by splitVoxels method
	 */
	static const unsigned short NO_SPLIT = 0;
	static const unsigned short NO_SPLIT_DUE_TO_VOXELS_LIMIT = 1;
	static const unsigned short NO_SPLIT_BUT_INITIALIZED = 2;
	static const unsigned short NORMAL_SPLIT = 3;

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


	virtual unsigned short splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc);

	// counts active top level voxels
	size_t countActiveTopLevelVoxels();

	virtual void getRandomEntry(RSAVector *position, RSAOrientation *orientation, Voxel **v, RND *rnd);
	Voxel *getVoxel(int i);
	virtual Voxel *getVoxel(const RSAVector &pos, const RSAOrientation &angle);
	virtual void getRandomPositionAndOrientation(RSAVector *position, RSAOrientation *orientation, Voxel *v, RND *rnd);
	double getSpatialVoxelSize() const;
	double getAngularVoxelSize() const;
	Voxel* get(int i);
	size_t getLength() const;
	virtual double getVoxelsVolume();
	std::string toPovray();
	std::string toWolfram();

	void store(std::ostream &f) const;
	void restore(std::istream &f);

};


#endif /* VOXELLIST_H_ */
