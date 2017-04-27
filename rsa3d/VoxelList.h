/*
 * VoxelList.h
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#ifndef VOXELLIST_H_
#define VOXELLIST_H_

#include "Voxel.h"
#include "RND.h"
#include "NeighbourGrid.h"
#include "BoundaryConditions.h"
#include "Positioned.h"
#include <unordered_set>

class VoxelList {

private:

	const double dxFactor = 1.001;
	NeighbourGrid* voxelNeighbourGrid;

	void fillNeighbourGrid();
	bool analyzeVoxel(Voxel *v, NeighbourGrid *nl, std::unordered_set<Positioned *> *neighbours, BoundaryConditions *bc);

	double *da;
	int* in;


protected:
	Voxel** voxels;
	int last;
	double voxelSize;
	double size;
	int beginningVoxelNumber;
	int dimension;

	Voxel* createVoxel(double* center, double vs, int index);
	void initVoxels(int n);
	void checkIndexes();

public:
	VoxelList(int N, double s, double d);

	virtual ~VoxelList();

	std::unordered_set<Positioned *> * getNeighbours(Voxel *v);
	void remove(Voxel *v);
	bool analyzeVoxel(Voxel *v, NeighbourGrid *nl, BoundaryConditions *bc, int timestamp);
	bool analyzeVoxel(Voxel *v, NeighbourGrid *nl, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel *v, std::unordered_set<Positioned *> *neighbours, BoundaryConditions *bc);
	bool splitVoxels(double minDx, int maxVoxels, NeighbourGrid *nl, BoundaryConditions *bc);

	Voxel *getRandomVoxel(RND *rnd);
	double* getRandomPosition(double *result, Voxel *v, RND *rnd);
	double getVoxelSize();
	Voxel* get(int i);
	int length();
	double getVoxelsSurface();




};

#endif /* VOXELLIST_H_ */
