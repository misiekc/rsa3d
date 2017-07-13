/*
 * OrientedCuboid.cpp
 *
 *  Created on: 12.07.2017
 *      Author: ciesla
 */

#include "OrientedCuboid.h"

bool OrientedCuboid::do2Drotation;

OrientedCuboid::OrientedCuboid(const Matrix & rotation) : Cuboid(rotation){
}

OrientedCuboid::~OrientedCuboid() {
}

void OrientedCuboid::initClass(const std::string &args){
    Cuboid::initClass(args);
    if (args.compare(args.length() - 4, 4, "true")==0)
    	OrientedCuboid::do2Drotation = true;
    else
    	OrientedCuboid::do2Drotation = false;
}

// Method creating (dynamically alocated) cuboid with random orientation.
// Used by ShapeFactory for shape generating
//----------------------------------------------------------------------------
Shape * OrientedCuboid::create(RND *rnd){
	OrientedCuboid *cuboid;
	if (OrientedCuboid::do2Drotation)
	    cuboid = new OrientedCuboid(Matrix::rotation3D(
	        0,
	        std::asin(rnd->nextValue() * 2 - 1),
	        rnd->nextValue() * 2 * M_PI));
	else
    	cuboid = new OrientedCuboid(Matrix::identity(3));

#ifdef CUBOID_DEBUG
    std::cout << "Creating OrientedCuboid:" << std::endl;
    std::cout << cuboid->orientation;
#endif

    return cuboid;
}

