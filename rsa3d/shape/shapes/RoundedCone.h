/*
 * RoundedCone.h
 *
 *  Created on: Aug 5, 2022
 *      Author: ciesla
 */

#ifndef SHAPE_SHAPES_ROUNDEDCONE_H_
#define SHAPE_SHAPES_ROUNDEDCONE_H_

#include "../ConvexShape.h"
#include "../../geometry/xenocollide/CollideGeometry.h"
#include "../../geometry/xenocollide/BodyBuilder.h"


class RoundedCone: public ConvexShape<3, 0> {

private:
	static MapPtr<CollideGeometry> shapeModel;
	static MapPtr<CollideGeometry> evModel;

	static double R, r, l;

	Matrix<3, 3> orientation;		// Cartesian rotation matrix 3D

public:
    static void initClass(const std::string &attr);
	explicit RoundedCone(const Matrix<3, 3> &orientation);
	virtual ~RoundedCone();

    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const override;
    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     double orientationRange) const override;
    Shape<3, 0> *clone() const override;

    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    std::string toPovray() const;
    std::string toWolfram() const;



};

#endif /* SHAPE_SHAPES_ROUNDEDCONE_H_ */
