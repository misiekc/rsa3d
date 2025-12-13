/*
 * GenericXCShape.h
 *
 *  Created on: Aug 7, 2022
 *      Author: ciesla
 */

#ifndef SHAPE_SHAPES_GENERICXCSHAPE_H_
#define SHAPE_SHAPES_GENERICXCSHAPE_H_

#include "../ConvexShape.h"
#include "../../geometry/xenocollide/CollideGeometry.h"
#include "../../geometry/xenocollide/MapPtr.h"

class GenericXCShape: public ConvexShape<3, 0> {

protected:

	static MapPtr<CollideGeometry> shapeModel;
	static MapPtr<CollideGeometry> evModel;

	Matrix<3, 3> orientation;		// Cartesian rotation matrix 3D


public:
	explicit GenericXCShape(const Matrix<3, 3> &orientation);
	virtual ~GenericXCShape();

    static void initClass(const std::string &attr);

    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const override;
    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     const Orientation<0> &orientationRange) const override;
    Shape<3, 0> *clone() const override;

    virtual void store(std::ostream &f) const override;
    virtual void restore(std::istream &f) override;
};

#endif /* SHAPE_SHAPES_GENERICXCSHAPE_H_ */
