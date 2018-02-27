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

#define ANGULAR_DIMENSION 0

template <unsigned short DIMENSION>
class VoxelList {

private:

	const double dxFactor = 1.0; // 1.0000000001;
	NeighbourGrid<Voxel<DIMENSION, ANGULAR_DIMENSION>>* voxelNeighbourGrid;
	bool* activeTopLevelVoxels;

	double findInitialVoxelSize(double d);
	int getLinearNumberOfVoxels(double vs);
	void fillNeighbourGrid();
//	bool analyzeVoxelOLD(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, std::unordered_set<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc);
//	bool analyzeVoxelNEW(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, std::unordered_set<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel<DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, std::vector<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc);
	bool disabled;

	std::uniform_real_distribution<double> *spatialDistribution;
	std::uniform_real_distribution<double> *angularDistribution;


protected:
	Voxel<DIMENSION, ANGULAR_DIMENSION>** voxels;
	int last;

	double initialVoxelSize;
	double voxelSize;
	double angularSize;

	double size;
	int beginningVoxelNumber;
	int offset[(1 << DIMENSION)][DIMENSION]; // matrix of d-dimensional offsets to 2^d voxel vertices


	Voxel<DIMENSION, ANGULAR_DIMENSION>* createVoxel(double* center, double *orientation, int index);
	void initVoxels();
	void checkIndexes();

public:

	VoxelList(double s, double d);
	void disable();

	virtual ~VoxelList();

	void getNeighbours(std::vector<Voxel<DIMENSION, ANGULAR_DIMENSION> *> *result, Voxel<DIMENSION, ANGULAR_DIMENSION> *v);
	void remove(Voxel<DIMENSION, ANGULAR_DIMENSION> *v);
	void removeTopLevelVoxel(Voxel<DIMENSION, ANGULAR_DIMENSION> *v);
	bool analyzeVoxel(Voxel<DIMENSION, ANGULAR_DIMENSION> *v, Shape<DIMENSION> *s, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel<DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, BoundaryConditions *bc, int timestamp);
	bool analyzeVoxel(Voxel<DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel<DIMENSION, ANGULAR_DIMENSION> *v, std::vector<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc);
	bool splitVoxels(double minDx, int maxVoxels, NeighbourGrid<Shape<DIMENSION>> *nl, BoundaryConditions *bc);

	Voxel<DIMENSION, ANGULAR_DIMENSION> *getRandomVoxel(RND *rnd);
	Voxel<DIMENSION, ANGULAR_DIMENSION> * getVoxel(double* da);
	double* getRandomPosition(double *result, Voxel<DIMENSION, ANGULAR_DIMENSION> *v, RND *rnd);
	double getVoxelSize();
	Voxel<DIMENSION, ANGULAR_DIMENSION>* get(int i);
	int length();
	double getVoxelsSurface();
	std::string toPovray();
};

#include "VoxelList.tpp"

#endif /* VOXELLIST_H_ */
