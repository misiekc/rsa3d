/*
 * ShapeFactory.cpp
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#include "ShapeFactory.h"

#include "shapes/Sphere.h"
#include "shapes/Cuboid.h"
#include "shapes/OrientedCuboid.h"
#include "shapes/SpheroCylinder2D.h"

Shape<RSA_DIMENSION>* (*ShapeFactory::createShape)(RND *rnd);

void ShapeFactory::initShapeClass(const std::string &sClass, const std::string &attr) {
	if (sClass == "Sphere") {
		Sphere<RSA_DIMENSION>::initClass(attr);
		ShapeFactory::createShape = Sphere<RSA_DIMENSION>::create;
	}
    else if (sClass == "OrientedCuboid") {
        OrientedCuboid<RSA_DIMENSION>::initClass(attr);
        ShapeFactory::createShape = OrientedCuboid<RSA_DIMENSION>::create;
    }
#if RSA_DIMENSION == 2
    else if (sClass == "SpheroCylinder2D")
    {
        SpheroCylinder2D::initClass(attr);
        ShapeFactory::createShape = SpheroCylinder2D::create;
    }
#elif RSA_DIMENSION == 3
	else if (sClass == "Cuboid") {
		Cuboid::initClass(attr);
		ShapeFactory::createShape = Cuboid::create3D;
	} else if (sClass == "Cuboid2D") {
		Cuboid::initClass(attr);
		ShapeFactory::createShape = Cuboid::create2D;
	}
#endif
    else {
        std::cerr << "Unknown shape: " << sClass << " or wrong RSA_DIMENTSION: " << RSA_DIMENSION << std::endl;
        exit(1);
    }

}
