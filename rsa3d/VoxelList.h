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

	const double dxFactor = 1.0; // 1.0000000001;
	NeighbourGrid<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>>* voxelNeighbourGrid;
	bool* activeTopLevelVoxels;

	double findInitialVoxelSize(double d);
	int getLinearNumberOfVoxels(double vs);
	void fillNeighbourGrid();
//	bool analyzeVoxelOLD(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, std::unordered_set<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc);
//	bool analyzeVoxelNEW(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, std::unordered_set<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *neighbours, BoundaryConditions *bc);
	bool disabled;

	std::uniform_real_distribution<double> *spatialDistribution;
	std::uniform_real_distribution<double> *angularDistribution;


protected:
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>** voxels;
	int last;

	double initialVoxelSize;
	double voxelSize;
	double angularSize;

	double size;
	int beginningVoxelNumber;
	int offset[(1 << SPATIAL_DIMENSION)][SPATIAL_DIMENSION]; // matrix of d-dimensional offsets to 2^d voxel vertices


	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>* createVoxel(double* center, double *orientation, int index);
	void initVoxels();
	void checkIndexes();

public:

	VoxelList(double s, double d);
	void disable();

	virtual ~VoxelList();

	void getNeighbours(std::vector<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *result, Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v);
	void remove(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v);
	void removeTopLevelVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v);
	bool analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc, int timestamp);
	bool analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *neighbours, BoundaryConditions *bc);
	bool splitVoxels(double minDx, int maxVoxels, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc);

	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *getRandomVoxel(RND *rnd);
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> * getVoxel(double* pos, double *angle);
	double* getRandomPositionAndOrientation(double *result, Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, RND *rnd);
	double getVoxelSize();
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>* get(int i);
	int length();
	double getVoxelsSurface();
	std::string toPovray();
	std::string toWolfram();
};

#include "VoxelList.tpp"

#endif /* VOXELLIST_H_ */
