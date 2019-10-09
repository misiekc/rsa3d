/*
 * OrientedCuboid.h
 *
 *  Created on: 12.07.2017
 *      Author: ciesla
 */

#ifndef SHAPES_ORIENTEDCUBOID_H_
#define SHAPES_ORIENTEDCUBOID_H_

#include "../../Voxel.h"
#include "../ConvexShape.h"


template <unsigned short DIMENSION>
class OrientedCuboid : public ConvexShape<DIMENSION, 0>{
private:
//	static bool do2Drotation;
    static double 			size[DIMENSION];

    int voxelFullyOrPartiallyInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &voxelPosition, double spatialSize) const;
protected:
    using EarlyRejectionResult = typename Shape<DIMENSION, 0>::EarlyRejectionResult;

public:
	static void initClass(const std::string &args);
    static bool voxelInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &voxelPosition,
    		double spatialSize, const std::vector<const RSAShape *> *shapes);

    bool overlap(BoundaryConditions<DIMENSION> *bc, const Shape<DIMENSION, 0> *s) const override;
	bool pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &position,
					 const Orientation<0> &orientation, double orientationRange) const override;
	bool voxelInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &voxelPosition,
					 const Orientation<0> &orientation, double spatialSize, double angularSize) const override;

    double getVolume(unsigned short dim) const override;
	Shape<DIMENSION, 0> *clone() const override;

	std::string toPovray() const;
};

#include "OrientedCuboid.tpp"

#endif /* SHAPES_ORIENTEDCUBOID_H_ */
