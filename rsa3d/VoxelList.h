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
#include <unordered_set>

template <ushort DIMENSION>
class VoxelList {

private:

	const double dxFactor = 1.0; // 1.0000000001;
	NeighbourGrid<Voxel<DIMENSION>>* voxelNeighbourGrid;
	bool* activeTopLevelVoxels;

	double findInitialVoxelSize(double d);
	int getLinearNumberOfVoxels(double vs);
	void fillNeighbourGrid();
//	bool analyzeVoxelOLD(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, std::unordered_set<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc);
//	bool analyzeVoxelNEW(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, std::unordered_set<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, std::unordered_set<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc);
	bool disabled;

	std::uniform_real_distribution<double> *distribution;


protected:
	Voxel<DIMENSION>** voxels;
	int last;

	double initialVoxelSize;
	double voxelSize;

	double size;
	int beginningVoxelNumber;
	int offset[(1 << DIMENSION)][DIMENSION]; // matrix of d-dimensional offsets to 2^d voxel vertices


	Voxel<DIMENSION>* createVoxel(double* center, double vs, int index);
	void initVoxels();
	void checkIndexes();

public:

	VoxelList(double s, double d);
	void disable();

	virtual ~VoxelList();

	void getNeighbours(std::unordered_set<Voxel<DIMENSION> *> *result, Voxel<DIMENSION> *v);
	void remove(Voxel<DIMENSION> *v);
	void removeTopLevelVoxel(Voxel<DIMENSION> *v);
	bool analyzeVoxel(Voxel<DIMENSION> *v, Shape<DIMENSION> *s, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, BoundaryConditions *bc, int timestamp);
	bool analyzeVoxel(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, BoundaryConditions *bc);
	bool analyzeVoxel(Voxel<DIMENSION> *v, std::unordered_set<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc);
	bool splitVoxels(double minDx, int maxVoxels, NeighbourGrid<Shape<DIMENSION>> *nl, BoundaryConditions *bc);

	Voxel<DIMENSION> *getRandomVoxel(RND *rnd);
	Voxel<DIMENSION> * getVoxel(double* da);
	double* getRandomPosition(double *result, Voxel<DIMENSION> *v, RND *rnd);
	double getVoxelSize();
	Voxel<DIMENSION>* get(int i);
	int length();
	double getVoxelsSurface();
	std::string toPovray();
};

#include "VoxelList.tpp"

#endif /* VOXELLIST_H_ */
