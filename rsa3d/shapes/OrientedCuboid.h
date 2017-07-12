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
public:
	OrientedCuboid(const Matrix & rotation);
	virtual ~OrientedCuboid();

	static void initClass(const std::string &args);
	static Shape * create(RND *rnd);
};

#endif /* SHAPES_ORIENTEDCUBOID_H_ */
