/*
 * VoxelList.h
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#ifndef VOXELLIST_H_
#define VOXELLIST_H_

#include <random>

#include "Voxel.h"
#include "RND.h"

class VoxelList {

private:

	static const double dxFactor = 1.001;
	NeighbourGrid* voxelNeighbourGrid;

	void fillNeighbourGrid();
	bool analyzeVoxel(Voxel *v, NeighbourGrid *nl, std::vector<Shape*> *neighbours, BoundaryConditions *bc);


protected:
	Voxel* voxels;
	int last;
	double voxelSize;
	double size;
	int beginningVoxelNumber;
	int dimension;

	Voxel* createVoxel(double* center, double vs, int index);
	void initVoxels(int n);
	void checkIndexes();

public:
	VoxelList(double s, double d);
	VoxelList(int N, double s, double d);

	virtual ~VoxelList();

	Voxel* getNeighbours(Voxel *v);
	void remove(Voxel *v);
	bool analyzeVoxel(Voxel *v, NeighbourGrid *nl, BoundaryConditions *bc, int timestamp);
	bool analyzeVoxel(Voxel *v, NeighbourGrid *nl, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel *v, std::vector<Shape*> *neighbours, BoundaryConditions *bc);
	bool splitVoxels(double minDx, int maxVoxels, NeighbourGrid *nl, BoundaryConditions *bc);

	Voxel *getRandomVoxel(RND *rnd);
	double* getRandomPosition(double *result, Voxel *v, RND *rnd);
	double getVoxelSize();
	Voxel* get(int i);
	int length();
	double getVoxelsSurface();




};

#endif /* VOXELLIST_H_ */
