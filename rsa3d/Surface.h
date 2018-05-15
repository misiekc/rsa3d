/*
 * Surface.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include "Parameters.h"
#include "Shape.h"
#include <vector>
#include <unordered_set>
#include "RND.h"
#include "Voxel.h"
#include "NeighbourGrid.h"
#include "VoxelList.h"

class Surface : public BoundaryConditions{

protected:

	NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> *list;

	double size;
	int dimension;

	void vectorFreeBC(double* v);
	void vectorPeriodicBC(double* v);

public:

	Surface(int dim, double s, double ndx, double vdx);
	virtual ~Surface();

	void clear();
	void add(Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s);
	Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>* check(Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s);
	void getNeighbours(std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *result, double *da);
	Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> * getClosestNeighbour(double *da, std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *neighbours);
	NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> * getNeighbourGrid();
	double distance2(const double *a1, const double *a2);

	virtual double * getTranslation(double *result, const double *p1, const double *p2) = 0;
	virtual void vector(double *v) = 0;
	virtual void checkPosition(double *da);
	virtual double getArea() = 0;

//	void drawShapes(Graphics g, double scale);
//	void drawShapes(Graphics g, double scale, double[] ta);
};

#endif /* SURFACE_H_ */
