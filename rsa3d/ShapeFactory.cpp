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

Shape* (*ShapeFactory::createShape)(RND *rnd);

void ShapeFactory::initShapeClass(const std::string &sClass, const std::string &attr){
	if (sClass.compare("Sphere")==0){
		Sphere::initClass(attr);
		ShapeFactory::createShape = Sphere::create;

	} else if (sClass.compare("Cuboid")==0){
		Cuboid::initClass(attr);
		ShapeFactory::createShape = Cuboid::create;

	} else if (sClass.compare("OrientedCuboid")==0){
		OrientedCuboid::initClass(attr);
		ShapeFactory::createShape = OrientedCuboid::create;

	}
}
