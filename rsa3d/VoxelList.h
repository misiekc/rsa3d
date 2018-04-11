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

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
class VoxelList {

private:

	// allows voxels to overlap - only for testing purposes and normally should be set to 1
	const double dxFactor = 1.0; // 1.0000000001;

	// neighbour grid structure for voxels. Needed to quickly find a voxel using its location
	NeighbourGrid<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>>* voxelNeighbourGrid;

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

	// checks consistency of indexes of root voxels
	void checkTopLevelVoxels();

	void compactVoxelArray(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **list, int &endIndex);

	// checks if a voxel can be removed
	bool analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *neighbours, BoundaryConditions *bc);


protected:
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>** voxels;
	int offset[(1 << SPATIAL_DIMENSION)][SPATIAL_DIMENSION]; // matrix of d-dimensional offsets to 2^d voxel vertices
	int last;

	double initialVoxelSize;
	double initialAngularVoxelSize;
	double voxelSize;
	double angularVoxelSize;


	double angularSize;
	double size;

	int beginningVoxelNumber;


	bool isVoxelInsidePacking(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v);
	bool isVoxelInsideExclusionZone(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, double spatialSize, double angularSize, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s, BoundaryConditions *bc);
	bool isVoxelInsideExclusionZone(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, double spatialSize, double angularSize, std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *shapes, BoundaryConditions *bc);
	bool isVoxelInsideExclusionZone(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, double spatialSize, double angularSize, std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *shapes, BoundaryConditions *bc, unsigned short maxDepth);

	void splitVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, double spatialSize, double angularSize, Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **vRes);


public:

	VoxelList(double s, double d, double ad);
	void disable();

	virtual ~VoxelList();

	void getNeighbours(std::vector<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *result, double *da);
	void remove(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v);
	void remove(std::vector<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *vVoxels);
	void removeTopLevelVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v);

	bool analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc, unsigned short depth);
	bool analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *neighbours, BoundaryConditions *bc);

	bool splitVoxels(double minDx, int maxVoxels, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc);

	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *getRandomVoxel(RND *rnd);
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *getVoxel(int i);
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> * getVoxel(double* pos, const double *angle);
	void getRandomPositionAndOrientation(double *position, double *orientation, Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, RND *rnd);
	double getVoxelSize();
	double getVoxelAngularSize();
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>* get(int i);
	int length() const;
	double getVoxelsSurface();
	std::string toPovray();
	std::string toWolfram();

	void store(std::ostream &f) const;
	void restore(std::istream &f);

};

#include "VoxelList.tpp"

#endif /* VOXELLIST_H_ */
