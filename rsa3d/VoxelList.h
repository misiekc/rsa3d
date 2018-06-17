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
#include "NeighbourGrid.h"
#include "BoundaryConditions.h"
#include "Positioned.h"
#include <vector>

class VoxelList {

private:

	// allows voxels to overlap - only for testing purposes and normally should be set to 1
	const double dxFactor = 1.0; // 1.0000000001;

	// neighbour grid structure for voxels. Needed to quickly find a voxel using its location
	NeighbourGrid<Voxel>* voxelNeighbourGrid;

	// probability distrubutions for drawing position inside a voxel
	std::uniform_real_distribution<double> *spatialDistribution;
	std::uniform_real_distribution<double> *angularDistribution;

	// array of top level voxel. If a top level voxel becomes inactive (due to shape placement, all its child voxels becomes obsolete
	bool* activeTopLevelVoxels;

	// disables voxel list for debug purposes - typically is false
	bool disabled;

	// returns initial size of a vovel. It is not grater than d and be an integer power of 2 (due to numerical issues)
	double findFloorSize(double d);

	// returns initial size of a vovel. It is not smaller than d and be an integer power of 2 (due to numerical issues)
	double findCeilSize(double d);

	// returns number of elements of size cellSize needed to cover a desired range
	size_t findArraySize(double range, double cellSize);

	// fills neigbour grid with voxels
	void fillNeighbourGrid();

	// for a given voxel returns index of its root
	int getIndexOfTopLevelVoxel(double *da);

	// initialize voxels - used inside a constructor
	void initVoxels();

	// checks voxels indexes consistency
//	void checkIndexes();

	// checks if a top level voxel for voxel v is active (if not v should be removed
	bool isTopLevelVoxelActive(Voxel *v);

	// checks consistency of indexes of root voxels
	void checkTopLevelVoxels();


	// aeeay of voxels will not have NULLs between pointers to objects - they can appear when splitting or analyze voxels in parallel
	void compactVoxelArray(Voxel **list, int &endIndex);

protected:
	Voxel** voxels;
	size_t length;

	double initialVoxelSize;
	double initialAngularVoxelSize;
	double spatialVoxelSize;
	double angularVoxelSize;


	double angularRange;
	double spatialRange;

	size_t beginningVoxelNumber;


	bool isVoxelInsidePacking(Voxel *v);
	bool isVoxelInsideExclusionZone(Voxel *v, double spatialSize, double angularSize, std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *shapes, BoundaryConditions *bc, unsigned short depth = 0);

	void splitVoxel(Voxel *v, double spatialSize, double angularSize, Voxel **vRes);


public:

	/**
	 * @brief Constructor
	 * @param packingSpatialSize packing size to be covered by voxels
	 * @param requestedSpatialVoxelSize suggested initial size of a voxel. Initial size of a allocated voxels will not be larger than the requested one.
	 * @param shapeAngularRange typpically 2*M_PI. Can be smaller for shapes with rotational symmetry. For example for squares it should be M_PI/2.0, and for ellipses, sherocylinders or rectangles M_PI.
	 * @param requestedAngularVoxelSize suggested initial angular size of a voxel. Initial angular size of allocaced voxels will not be larger than requested one.
	 */
	VoxelList(double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange, double requestedAngularVoxelSize);

	void disable();

	virtual ~VoxelList();

	void getNeighbours(std::vector<Voxel *> *result, double *da);
	void remove(Voxel *v);
	void remove(std::vector<Voxel *> *vVoxels);
	void removeTopLevelVoxel(Voxel *v);

	bool analyzeVoxel(Voxel *v, Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s, BoundaryConditions *bc);
//	bool analyzeVoxel(Voxel *v, std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *neighbours, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel *v, NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc, unsigned short depth = 0);

	size_t analyzeVoxels(BoundaryConditions *bc, NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> *nl, unsigned short depth);


	bool splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc);

	Voxel *getRandomVoxel(RND *rnd);
	Voxel *getVoxel(int i);
	Voxel * getVoxel(double* pos, const std::array<double, RSA_ANGULAR_DIMENSION> &angle);
	void getRandomPositionAndOrientation(double *position, double *orientation, Voxel *v, RND *rnd);
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
