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

	// for a given voxel returns index of its root
	int getIndexOfTopLevelVoxel(const RSAVector &da) const;

	// sets status basing on existing voxels
	void refreshTopLevelVoxels();

	// finds voxel containing given point - not used - only for debugging
	Voxel * findVoxel(Voxel **list, size_t listSize, const RSAVector &pos,
					  const RSAOrientation &angle);
protected:
    // allows voxels to overlap - only for testing purposes and normally should be set to 1
    const double dxFactor = 1.0; // 1.0000000001;
    // disables voxel list for debug purposes - typically is false
	bool disabled;
	// to optimize memory occupations voxels are initialized in the first invoke split()
	bool voxelsInitialized;

	unsigned short int surfaceDimension;
	Voxel** voxels;
	size_t length;

	double initialVoxelSize;
	RSAOrientation initialAngularVoxelSize;

	double spatialVoxelSize;
	RSAOrientation angularVoxelSize;

	RSAOrientation angularRange;
	double spatialRange;

	// array of top level voxel. If a top level voxel becomes inactive (due to shape placement, all its child voxels becomes obsolete
	bool* activeTopLevelVoxels;

	// probability distrubutions for drawing position inside a voxel
	std::uniform_real_distribution<double> *spatialDistribution;
	#if RSA_ANGULAR_DIMENSION > 0
		std::uniform_real_distribution<double> *angularDistribution[RSA_ANGULAR_DIMENSION];
	#else
		std::uniform_real_distribution<double> **angularDistribution;
	#endif
	//size_t beginningVoxelNumber;

	// neighbour grid structure for voxels. Needed to quickly find a voxel using its location
	NeighbourGrid<Voxel>* voxelNeighbourGrid;

	VoxelList();

	// returns initial size of a vovel. It is not grater than d and be an integer power of 2 (due to numerical issues)
	static double findFloorSize(double d);

	// returns initial size of a vovel. It is not smaller than d and be an integer power of 2 (due to numerical issues)
	static double findCeilSize(double d);

	virtual void allocateVoxels(size_t size);

	// initialize voxels, returns number of initial voxels
	virtual unsigned int initVoxels(RSABoundaryConditions *bc, NeighbourGrid<const RSAShape> *nl);

	// checks consistency of indexes of root voxels
	void checkTopLevelVoxels();

	// checks if a top level voxel for voxel v is active (if not v should be removed
	bool isTopLevelVoxelActive(Voxel *v) const;

	virtual bool isVoxelInsidePacking(const Voxel *v, double spatialSize) const;
	virtual bool isVoxelInsideExclusionZone(Voxel *v, double spatialSize, RSAOrientation angularSize,
                                    std::vector<const RSAShape*> *shapes, RSABoundaryConditions *bc,
                                    unsigned short depth = 0) const;

	bool isVoxelInsideExclusionZoneOld(Voxel *v, double spatialSize, const RSAOrientation &angularSize,
                                    std::vector<const RSAShape*> *shapes, RSABoundaryConditions *bc,
                                    unsigned short depth = 0);

	virtual void splitVoxel(Voxel *v, double spatialSize, const RSAOrientation &angularSize, Voxel **vRes) const;

	// fills neigbour grid with voxels
	void rebuildNeighbourGrid();

	virtual bool analyzeVoxel(Voxel *v, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc, double spatialSize, RSAOrientation angularSize, unsigned short depth = 0) const;

	virtual void moveVoxelInList(size_t from, size_t to);

	// voxels array will not have NULLs between pointers to objects - they can appear when splitting or analyze voxels in parallel
	// returns number of not null values in list
	void compactVoxelArray();

	virtual void getRandomPositionAndOrientation(RSAVector *position, RSAOrientation *orientation, Voxel *v, RND *rnd);

	// Returns number of grid cells in each direction of the simulation box.
	// The default implementation gives the same number of cells in each directions covering the whole simulation box,
	// but derived classes can alter this behaviour
    [[nodiscard]] virtual std::array<std::size_t, RSA_SPATIAL_DIMENSION> calculateSpatialGridLinearSize() const;

    // returns number of elements of size cellSize needed to cover a desired range
    std::size_t findArraySize(double range, double cellSize) const;

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
	VoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, RSAOrientation shapeAngularRange, RSAOrientation requestedAngularVoxelSize);

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
	virtual Voxel *getVoxel(const RSAVector &pos, const RSAOrientation &angle) const;
	virtual double getSpatialVoxelSize() const;
	virtual RSAOrientation getAngularVoxelSize() const;
	Voxel* get(int i);
	size_t getLength() const;
	virtual double getVoxelsVolume() const;
	std::string toPovray() const;
	std::string toWolfram() const;

	void store(std::ostream &f) const;
	void restore(std::istream &f);
};


#endif /* VOXELLIST_H_ */
