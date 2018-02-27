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

template <unsigned short DIMENSION>
class Shape : public Positioned<DIMENSION>{

public:
	static Shape<DIMENSION>* (*createShape)(RND *rnd);

	int no;
	double time;

	Shape();
	Shape(const Shape<DIMENSION> & other);

	virtual ~Shape();

	// returns linear size of a cell in a NeighbourGrig. This size should be as small as possible but big enough to avoid overlapping between shapes having centers in cells that are not neighbours
	virtual double getNeighbourListCellSize() = 0;

	// returns initial linear size of a (cubic) voxel. This size should be as big as possible but shape with the center inside the voxel have to cover the whole voxel
	virtual double getVoxelSize() = 0;

	// translates the shape by the given vector v
	virtual void translate(double* v);

	// checks if there is overlap with the shape pointed by s.
	virtual int overlap(BoundaryConditions *bc, Shape *s) = 0;

	// returns a volume of the shape
	virtual double getVolume() = 0;

	// checks if a given point is within excluded volume
	virtual int pointInside(BoundaryConditions *bc, double* da) = 0;

	virtual int pointInside(BoundaryConditions *bc, double* position, double *orientation, double orientationRange);

	// moves the shape towards given shape s
	virtual double minDistance(Shape<DIMENSION> *s);

	// returns string representation of the shape
	virtual std::string toString();

	// returns povray string representation of the shape
	virtual std::string toPovray() const;

	// serialize shape
	virtual void store(std::ostream &f) const;

	// deserialize shape
	virtual void restore(std::istream &f);

	// draws the shape
    // virtual void draw()

    void vectorTranslate(const Vector<DIMENSION> &translation);
	void applyBC(BoundaryConditions *bc, Shape<DIMENSION> *second);
	Vector<2> applyBC(BoundaryConditions *bc, double *pointToTranslate);
};

#include "Shape.tpp"


#endif /* SHAPE_H_ */
