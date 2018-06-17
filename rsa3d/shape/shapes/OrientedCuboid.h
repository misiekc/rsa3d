/*
 * OrientedCuboid.h
 *
 *  Created on: 12.07.2017
 *      Author: ciesla
 */

#ifndef SHAPES_ORIENTEDCUBOID_H_
#define SHAPES_ORIENTEDCUBOID_H_


template <unsigned short DIMENSION>
class OrientedCuboid : public ConvexShape<DIMENSION, 0>{
private:
//	static bool do2Drotation;
    static double 			size[DIMENSION];

public:
	OrientedCuboid();

	virtual ~OrientedCuboid();
	static void initClass(const std::string &args);

    bool overlap(BoundaryConditions *bc, const Shape<DIMENSION, 0> *s) const override;
    bool pointInside(BoundaryConditions *bc, double* da) const override;
	bool pointInside(BoundaryConditions *bc, double* position, const std::array<double, 0> &orientation,
					double orientationRange) const override;
    double getVolume() const;
	Shape<DIMENSION, 0> *clone() const override;

	std::string toPovray() const;
};

#include "OrientedCuboid.tpp"

#endif /* SHAPES_ORIENTEDCUBOID_H_ */
