/*
 * OrientedCuboid.h
 *
 *  Created on: 12.07.2017
 *      Author: ciesla
 */

#ifndef SHAPES_ORIENTEDCUBOID_H_
#define SHAPES_ORIENTEDCUBOID_H_

#include "Cuboid.h"

class OrientedCuboid: public Cuboid {
private:
	static bool do2Drotation;

public:
	OrientedCuboid(const Matrix<3, 3> & rotation);
	virtual ~OrientedCuboid();

	static void initClass(const std::string &args);
	static Shape * create(RND *rnd);

    int overlap(BoundaryConditions *bc, Shape *s);
    int pointInside(BoundaryConditions *bc, double* da);
};

#endif /* SHAPES_ORIENTEDCUBOID_H_ */
