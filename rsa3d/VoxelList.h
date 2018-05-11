/*
 * VoxelList.h
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#ifndef VOXELLIST_H_
#define VOXELLIST_H_

#include "Voxel.h"
#include "Shape.h"
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

	// returns initial spatial size of a vovel. It should be not grater than d and be an integer power of 2 (due to numerical issues)
	double findInitialVoxelSize(double d);

	// returns initial angular size of a vovel. It should be not smaller than d and be an integer power of 2 (due to numerical issues)
	double findInitialVoxelAngularSize(double d);

	// returns number ov voxels needed to cover a packing (along a line)
	int getLinearNumberOfVoxels(double vs);

	// fills neigbour grid with voxels
	void fillNeighbourGrid();

	// for a given voxel returns index of its root
	int getIndexOfTopLevelVoxel(double *da);

	// initialize voxels - used inside a constructor
	void initVoxels();

	// checks voxels indexes consistency
	void checkIndexes();

	// checks if a top level voxel for voxel v is active (if not v should be removed
	bool isTopLevelVoxelActive(Voxel *v);

	// checks consistency of indexes of root voxels
	void checkTopLevelVoxels();

	// aeeay of voxels will not have NULLs between pointers to objects - they can appear when splitting or analyze voxels in parallel
	void compactVoxelArray(Voxel **list, int &endIndex);

protected:
	Voxel** voxels;
	int offset[(1 << RSA_SPATIAL_DIMENSION)][RSA_SPATIAL_DIMENSION]; // matrix of d-dimensional offsets to 2^d voxel vertices
	int last;

	double initialVoxelSize;
	double initialAngularVoxelSize;
	double voxelSize;
	double angularVoxelSize;


	double angularSize;
	double size;

	int beginningVoxelNumber;


	bool isVoxelInsidePacking(Voxel *v);
	bool isVoxelInsideExclusionZone(Voxel *v, double spatialSize, double angularSize, Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s, BoundaryConditions *bc);
	bool isVoxelInsideExclusionZone(Voxel *v, double spatialSize, double angularSize, std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *shapes, BoundaryConditions *bc, unsigned short depth = 0);

	void splitVoxel(Voxel *v, double spatialSize, double angularSize, Voxel **vRes);


public:

	VoxelList(double s, double d, double ad);
	void disable();

	virtual ~VoxelList();

	void getNeighbours(std::vector<Voxel *> *result, double *da);
	void remove(Voxel *v);
	void remove(std::vector<Voxel *> *vVoxels);
	void removeTopLevelVoxel(Voxel *v);

	bool analyzeVoxel(Voxel *v, Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s, BoundaryConditions *bc);
//	bool analyzeVoxel(Voxel *v, std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *neighbours, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel *v, NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc, unsigned short depth = 0);

	bool splitVoxels(double minDx, int maxVoxels, NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc);

	Voxel *getRandomVoxel(RND *rnd);
	Voxel *getVoxel(int i);
	Voxel * getVoxel(double* pos, const std::array<double, RSA_ANGULAR_DIMENSION> &angle);
	void getRandomPositionAndOrientation(double *position, double *orientation, Voxel *v, RND *rnd);
	double getVoxelSize();
	double getVoxelAngularSize();
	Voxel* get(int i);
	int length() const;
	double getVoxelsSurface();
	std::string toPovray();
	std::string toWolfram();

	void store(std::ostream &f) const;
	void restore(std::istream &f);

};


#endif /* VOXELLIST_H_ */
