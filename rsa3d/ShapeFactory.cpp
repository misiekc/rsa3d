/*
 * ShapeFactory.cpp
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#include "ShapeFactory.h"

#include "shapes/Sphere.h"

Shape* (*ShapeFactory::createShape)(RND *rnd);

void ShapeFactory::initShapeClass(const std::string &sClass, const std::string &attr){
	if (sClass.compare("Sphere")==0){
		Sphere::initClass(attr);
		ShapeFactory::createShape = Sphere::create;

	}
}
