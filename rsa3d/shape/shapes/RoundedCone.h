/*
 * RoundedCone.h
 *
 *  Created on: Aug 5, 2022
 *      Author: ciesla
 */

#ifndef SHAPE_SHAPES_ROUNDEDCONE_H_
#define SHAPE_SHAPES_ROUNDEDCONE_H_

#include "GenericXCShape.h"


class RoundedCone: public GenericXCShape {

private:
	static double R, r, l;
	static double volume();


public:
    static void initClass(const std::string &attr);
	explicit RoundedCone(const Matrix<3, 3> &orientation);
	virtual ~RoundedCone();

//    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const override;
//    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation, double orientationRange) const override;
    Shape<3, 0> *clone() const override;

    std::string toPovray() const;
    std::string toWolfram() const;
};

#endif /* SHAPE_SHAPES_ROUNDEDCONE_H_ */
