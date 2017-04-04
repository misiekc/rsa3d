/*
 * Surface.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include "Shape.h"
#include <vector>
#include "RND.h"

class Surface : public BoundaryConditions{

private:
	static const int FACTOR_LIMIT = 5;

	std::vector<Shape *> shapes;
	int missCounter;
	int seed=-1;

	int tmpSplit, iAnalyze, iMaxVoxels;
	double dMinVoxelSize;


	int analyzeVoxels();
	int analyzeRegion(Voxel *v);


protected:

	NeighbourGrid *list;
	// VoxelList *voxels;
	double size;
	int dimension;

	void vectorFreeBC(double* v);
	void vectorPeriodicBC(double* v);




public:
	Surface(int dim, double s, double ndx, double vdx);
	virtual ~Surface();

	void setSeed(int s);
	void add(Shape *s);
	bool check(Shape *s);
	std::vector<Shape *> getNeighbours(double *da);
	double distance2(double *a1, double *a2);

	virtual double * getTranslation(double *result, double *p1, double *p2) = 0;
	virtual void vector(double *v) = 0;
	virtual double * getRandomPosition(double *result, RND *rnd) = 0;
	virtual double getArea() = 0;

	void setParameters(int ia, int is, double dvs, int imv);
	bool doIteration(Shape *s, RND *rnd);

	bool isSaturated();

	double getFactor();

	std::vector<Shape *> getShapes();

//	void drawShapes(Graphics g, double scale);
//	void drawShapes(Graphics g, double scale, double[] ta);
};

#endif /* SURFACE_H_ */
