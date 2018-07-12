/*
 * ShapeFactory.cpp
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#include "ShapeFactory.h"
#include "../Positioned.h"
#include "../Utils.h"

#include "shapes/Sphere.h"
#include "shapes/cuboid/Cuboid.h"
#include "shapes/OrientedCuboid.h"
#include "shapes/SpheroCylinder2D.h"
#include "shapes/Ellipse.h"
#include "shapes/Rectangle.h"
#include "shapes/polygon/Polygon.h"
#include "shapes/polygon/SBPolygon.h"
#include "shapes/polygon/HBPolygon.h"
#include "shapes/regular_solid/Dodecahedron.h"
#include "shapes/regular_solid/Icosahedron.h"
#include "shapes/regular_solid/Octahedron.h"
#include "shapes/regular_solid/Tetrahedron.h"
#include "shapes/regular_solid/Cuboctahedron.h"
#include "shapes/regular_solid/TruncatedCube.h"
#include "shapes/regular_solid/TruncatedOctahedron.h"
#include "shapes/regular_solid/Rhombicuboctahedron.h"

void ShapeFactory::initShapeClass(const std::string &sClass, const std::string &attr) {

	int in[RSA_SPATIAL_DIMENSION];
	for(ushort i=0; i<RSA_SPATIAL_DIMENSION; i++)
		in[i] = 0;
	int index = 0;
	do{
		std::copy(in, in+RSA_SPATIAL_DIMENSION, Positioned<RSA_SPATIAL_DIMENSION>::offset[index]);
		index++;
	}while(increment(in, RSA_SPATIAL_DIMENSION, 1));


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
            } else if (sClass == "Rectangle") {
                Rectangle::initClass(attr);
                return;
            } else if (sClass == "SBPolygon") {
            	SBPolygon::initClass(attr);
            	return;
            } else if (sClass == "HBPolygon") {
            	HBPolygon::initClass(attr);
            	return;
            } else if (sClass == "Polygon") {
            	Polygon::initClass(attr);
            	return;
            }
        #endif

    #elif RSA_SPATIAL_DIMENSION == 3 && RSA_ANGULAR_DIMENSION == 0
        if (sClass == "Cuboid") {
            Cuboid::initClass(attr);
            return;
        } else if (sClass == "Tetrahedron") {
            RegularSolid<Tetrahedron>::initClass(attr);
            return;
        } else if (sClass == "Octahedron") {
            RegularSolid<Octahedron>::initClass(attr);
            return;
        } else if (sClass == "Dodecahedron") {
            RegularSolid<Dodecahedron>::initClass(attr);
            return;
        } else if (sClass == "Icosahedron") {
            RegularSolid<Icosahedron>::initClass(attr);
            return;
        } else if (sClass == "Cuboctahedron") {
            RegularSolid<Cuboctahedron>::initClass(attr);
            return;
        } else if (sClass == "TruncatedCube") {
            RegularSolid<TruncatedCube>::initClass(attr);
            return;
        } else if (sClass == "TruncatedOctahedron") {
            RegularSolid<TruncatedOctahedron>::initClass(attr);
            return;
        } else if (sClass == "Rhombicuboctahedron") {
            RegularSolid<Rhombicuboctahedron>::initClass(attr);
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
