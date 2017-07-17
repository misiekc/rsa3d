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
#include <unordered_set>
#include "RND.h"
#include "Voxel.h"
#include "NeighbourGrid.h"
#include "VoxelList.h"

class Surface : public BoundaryConditions{

protected:

	NeighbourGrid *list;

	double size;
	int dimension;

	void vectorFreeBC(double* v);
	void vectorPeriodicBC(double* v);

public:

	Surface(int dim, double s, double ndx, double vdx);
	virtual ~Surface();

	void add(Shape *s);
	Shape* check(Shape *s);
	std::unordered_set<Positioned *> * getNeighbours(double *da);
	NeighbourGrid * getNeighbourGrid();
	double distance2(double *a1, double *a2);

	virtual double * getTranslation(double *result, double *p1, double *p2) = 0;
	virtual void vector(double *v) = 0;
	virtual double * getRandomPosition(double *result, RND *rnd) = 0;
	virtual double getArea() = 0;

//	void drawShapes(Graphics g, double scale);
//	void drawShapes(Graphics g, double scale, double[] ta);
};

#endif /* SURFACE_H_ */
