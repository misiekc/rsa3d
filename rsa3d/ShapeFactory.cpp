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
#include "shapes/Ellipse.h"
#include "shapes/Rectangle.h"

void ShapeFactory::initShapeClass(const std::string &sClass, const std::string &attr) {

    // Shapes for any dimension, without angular dimensions
    #if RSA_ANGULAR_DIMENSION == 0
        if (sClass == "Sphere") {
            Sphere<RSA_SPATIAL_DIMENSION>::initClass(attr);
            return;
        } else if (sClass == "OrientedCuboid") {
            OrientedCuboid<RSA_SPATIAL_DIMENSION>::initClass(attr);
            return;
        }
    #endif

    // Shapes of specific dimensions
    #if RSA_SPATIAL_DIMENSION == 2
        // 2D shapes with angular dimension
        #if RSA_ANGULAR_DIMENSION == 1
            if (sClass == "SpheroCylinder2D") {
                SpheroCylinder2D::initClass(attr);
                return;
            } else if (sClass == "Ellipse") {
                Ellipse::initClass(attr);
                return;
 //           } else if (sClass == "Rectangle") {
 //               Rectangle::initClass(attr);
 //               ShapeFactory::createShape = Rectangle::create;
 //               return;
            }
        #endif

    #elif RSA_SPATIAL_DIMENSION == 3 && RSA_ANGULAR_DIMENSION == 0
        if (sClass == "Cuboid") {
            Cuboid::initClass(attr);
            return;
        }
    #endif

    std::cerr << "Unknown shape: " << sClass << " or wrong dimensions: " << RSA_SPATIAL_DIMENSION << ", " << RSA_ANGULAR_DIMENSION << std::endl;
    exit(EXIT_FAILURE);
}

RSAShape *ShapeFactory::createShape(RND *rnd) {
    // Fetch appropriate function from Shape
    auto createShapeImpl = RSAShape::getCreateShapeImpl();
    return createShapeImpl(rnd);
}
