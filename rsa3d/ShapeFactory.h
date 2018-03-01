/*
 * ShapeFactory.h
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#ifndef SHAPEFACTORY_H_
#define SHAPEFACTORY_H_

#include "Parameters.h"
#include "Shape.h"
#include <string>

class ShapeFactory {
public:
	static Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>* (*createShape)(RND *rnd);
	static void initShapeClass(const std::string &sClass, const std::string &attr);
};

#endif /* SHAPEFACTORY_H_ */
