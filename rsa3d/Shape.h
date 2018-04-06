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
#include <array>

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
class Shape : public Positioned<SPATIAL_DIMENSION>{

private:

    static double voxelAngularSize;

    std::array<double, ANGULAR_DIMENSION> orientation;

protected:

    /**
     * Translates second by a vector returned by
     * \code
     * bc->getTranslation(result, this, second);
     * \endcode
     * @param bc boundary conditions to apply
     * @param second shape to be translated according to bc
     */
	virtual void applyBC(BoundaryConditions *bc, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *second) const;

    /**
     * Sets shape's orientation. Derived shape classes have to override this method when they want to keep track of
     * shape's orientation. One would then typically write:
     * \code
     * void Derived::setOrientation(const double *orientation) {
     *     Shape::setOrientation(orientation);
     *     // invoke getOrientation and for example compute vertices
     * }
     * \endcode
     * rotate(double*) method delegates to this, so no orientation changes will be missed.
     * @param orientation new shape's orientation
     */
    virtual void setOrientation(const double *orientation);

public:

	static Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>* (*createShape)(RND *rnd);
    int no;

	double time;

    Shape();
	Shape(const Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> & other);
    Shape &operator=(const Shape &other);
	virtual ~Shape();

    /**
     * Returns linear size of a cell in a NeighbourGrid. This size should be as small as possible but big enough to
     * avoid overlapping between shapes having centers in cells that are not neighbours.
     * @return linear size of a cell in a NeighbourGrid
     */
	virtual double getNeighbourListCellSize() const = 0;

    /**
     * Returns initial linear size of a (cubic) voxel. This size should be as big as possible but shape with the center
     * inside the voxel have to cover the whole voxel.
     * @return initial linear size of a (cubic) voxel
     */
	virtual double getVoxelSize() const = 0;

    /**
     * Returns angular size of a voxel. This should be as small as possible but big enough so that angle range
     * (0, `getVoxelAngularSize()`) describes all possible shape's orientations.
     * @return angular size of a voxel
     */
	virtual double getVoxelAngularSize() const;

    /**
     * Returns an array of all angles describing shape's orientation.
     * @return array describing shape's orietation
     */
    const double * getOrientation() const;

	/**
	 * Translates the shape by a given vector v.
	 * @param v a vector to translate
	 */
    void translate(double* v);

    /**
     * Increases all shape's angles by respective values from array v.
     * @param v an array of angle deltas
     */
    void rotate(double *v);

	/**
	 * Checks if there is overlap with the shape pointed by s.
	 * @param bc boundary conditions to take into account
	 * @param s the second shape
	 * @return 0 if there is no overlap, nonzero number otherwise
	 */
	virtual int overlap(BoundaryConditions *bc, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s) const = 0;

    /**
     * Returns a volume of the shape.
     * @return a volume of the shape
     */
	virtual double getVolume() const = 0;

    /**
     * Checks if a given point is within intersection of all excluded volumes for all orientations denoted by
     * orientation and orientationRange parameters.
     * @param bc boundary conditions to take into account
     * @param position position of a point to check
     * @param orientation array of beginnings of angle intervals
     * @param orientationRange array of lengths of angle intervals
     * @return 0 if point is outside, nonzero number otherwise
     */
	virtual int pointInside(BoundaryConditions *bc, double* position, double *orientation, double orientationRange) const = 0;

    /**
     * Checks if a given point is within excluded volume for any orientation of a virtual rotationally symmetric shape.
     * @param bc boundary conditions to take into account
     * @param da position of a point to check
     * @return 0 if point is outside, nonzero number otherwise
     */
	virtual int pointInside(BoundaryConditions *bc, double* da) const;

    /**
     * ?? Moves the shape towards given shape s.
     * @param s ??
     * @return ??
     */
	virtual double minDistance(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s) const;

    /**
     * Returns string representation of the shape.
     * @return string representation of the shape
     */
	virtual std::string toString() const;

    /**
     * Returns povray string representation of the shape.
     * @return povray string representation of the shape
     */
	virtual std::string toPovray() const;

    /**
     * Returns Wolfram Mathematica string representation of the shape.
     * @return Wolfram Mathematica string representation of the shape
     */
	virtual std::string toWolfram() const;

    /**
     * Serializes shape and writes it in binary form to f.
     * @param f stream to write serialized shape to
     */
	virtual void store(std::ostream &f) const;

	/**
	 * Deserializes shape from f.
	 * @param f stream to write serialized shape from
	 */
	virtual void restore(std::istream &f);
};

#include "Shape.tpp"


#endif /* SHAPE_H_ */
