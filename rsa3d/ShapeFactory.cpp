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

Shape<RSA_DIMENSION>* (*ShapeFactory::createShape)(RND *rnd);

void ShapeFactory::initShapeClass(const std::string &sClass, const std::string &attr){
	if (sClass.compare("Sphere")==0){
		Sphere<RSA_DIMENSION>::initClass(attr);
		ShapeFactory::createShape = Sphere<RSA_DIMENSION>::create;

	}

#if RSA_DIMENSION == 3
	else if (sClass.compare("Cuboid")==0){
		Cuboid::initClass(attr);
		ShapeFactory::createShape = Cuboid::create;

	}
#endif

	else if (sClass.compare("OrientedCuboid")==0){
		OrientedCuboid<RSA_DIMENSION>::initClass(attr);
		ShapeFactory::createShape = OrientedCuboid<RSA_DIMENSION>::create;

	}
}
