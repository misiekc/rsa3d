/*
 * OrientedCuboid.cpp
 *
 *  Created on: 12.07.2017
 *      Author: ciesla
 */

#include "OrientedCuboid.h"

OrientedCuboid::OrientedCuboid(const Matrix & rotation) : Cuboid(rotation){
}

OrientedCuboid::~OrientedCuboid() {
}

void OrientedCuboid::initClass(const std::string &args)
{
    Cuboid::initClass(args);
}

// Method creating (dynamically alocated) cuboid with random orientation.
// Used by ShapeFactory for shape generating
//----------------------------------------------------------------------------
Shape * OrientedCuboid::create(RND *rnd)
{
    OrientedCuboid * cuboid = new OrientedCuboid(Matrix::identity(3));

#ifdef CUBOID_DEBUG
    std::cout << "Creating OrientedCuboid:" << std::endl;
    std::cout << cuboid->orientation;
#endif

    return cuboid;
}

