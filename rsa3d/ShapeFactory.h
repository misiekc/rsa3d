/*
 * ShapeFactory.h
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#ifndef SHAPEFACTORY_H_
#define SHAPEFACTORY_H_

#include "Shape.h"
#include <string>

class ShapeFactory {
public:
	static Shape* (*createShape)(RND *rnd);
	static void initShapeClass(const std::string &sClass, const std::string &attr);
};

#endif /* SHAPEFACTORY_H_ */
