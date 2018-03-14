/*
 * Shape.h
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#ifndef SHAPE_H_
#define SHAPE_H_

#include "BoundaryConditions.h"
#include "Positioned.h"
#include "RND.h"
#include <string>
#include <ostream>
#include <istream>

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
class Shape : public Positioned<SPATIAL_DIMENSION>{

private:
	static double voxelAngularSize;

protected:
	double orientation[ANGULAR_DIMENSION];

	void applyBC(BoundaryConditions *bc, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *second);

public:

	static Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>* (*createShape)(RND *rnd);

	int no;
	double time;

	Shape();
	Shape(const Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> & other);

	virtual ~Shape();

	// returns linear size of a cell in a NeighbourGrig. This size should be as small as possible but big enough to avoid overlapping between shapes having centers in cells that are not neighbours
	virtual double getNeighbourListCellSize() = 0;

	// returns initial linear size of a (cubic) voxel. This size should be as big as possible but shape with the center inside the voxel have to cover the whole voxel
	virtual double getVoxelSize() = 0;

	virtual double getVoxelAngularSize();

	double * getOrientation();

	// translates the shape by the given vector v
	void translate(double* v);

	virtual void rotate(double *v);

	// checks if there is overlap with the shape pointed by s.
	virtual int overlap(BoundaryConditions *bc, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s) = 0;

	// returns a volume of the shape
	virtual double getVolume() = 0;

	// checks if a given point is within excluded volume
	virtual int pointInside(BoundaryConditions *bc, double* position, double *orientation, double orientationRange) = 0;

	// checks if a given point is within excluded volume for any orientation of a virtual rotationally symmetric shape
	virtual int pointInside(BoundaryConditions *bc, double* da);

	// moves the shape towards given shape s
	virtual double minDistance(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s);

	// returns string representation of the shape
	virtual std::string toString();

	// returns povray string representation of the shape
	virtual std::string toPovray() const;

	// returns Wolfram Mathematica string representation of the shape
	virtual std::string toWolfram() const;

	// serialize shape
	virtual void store(std::ostream &f) const;

	// deserialize shape
	virtual void restore(std::istream &f);

	// draws the shape
    // virtual void draw()

//    void vectorTranslate(const Vector<SPATIAL_DIMENSION, ANGULAR_DIMENSION> &translation);

};

#include "Shape.tpp"


#endif /* SHAPE_H_ */
