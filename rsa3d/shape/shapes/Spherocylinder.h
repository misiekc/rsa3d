/*
 * Spherocylinder.h
 *
 *  Created on: Mar 15, 2019
 *      Author: ciesla
 */

#ifndef SHAPE_SHAPES_SPHEROCYLINDER_H_
#define SHAPE_SHAPES_SPHEROCYLINDER_H_

#include "../ConvexShape.h"
#include "../OrderCalculable.h"


template <unsigned short DIMENSION>
class Spherocylinder : public ConvexShape<DIMENSION, 0> , public OrderCalculable {

private:
	static double lenght;
	static double radius;

	Matrix<DIMENSION, DIMENSION> orientation;

	void getEnd(short beginOrEnd, Vector<DIMENSION> *result) const;

	double distanceFrom(Spherocylinder *s) const;
	double distanceFrom(Vector<DIMENSION> *v) const;


public:
    static void initClass(const std::string &attr);

    explicit Spherocylinder(const Matrix<DIMENSION, DIMENSION> &orientation);

    bool overlap(BoundaryConditions<DIMENSION> *bc, const Shape<DIMENSION, 0> *s) const override;
    bool pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &position, const Orientation<0> &orientation, double orientationRange) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    Shape<DIMENSION, 0> *clone() const override;

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;

    std::string toPovray() const;
//    std::string toWolfram() const;
};

#include "Spherocylinder.tpp"

#endif /* SHAPE_SHAPES_SPHEROCYLINDER_H_ */
