/*
 * ShapeFactory.cpp
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#include "ShapeFactory.h"
#include "../Positioned.h"
#include "../utils/Assertions.h"

#if RSA_ANGULAR_DIMENSION == 0
	#include "shapes/Sphere.h"
	#include "shapes/spherocylinder/Spherocylinder.h"
	#include "shapes/OrientedCuboid.h"

	#if RSA_SPATIAL_DIMENSION == 3
		#include "shapes/cuboid/Cuboid.h"
		#include "shapes/regular_solid/Dodecahedron.h"
		#include "shapes/regular_solid/Icosahedron.h"
		#include "shapes/regular_solid/Octahedron.h"
		#include "shapes/regular_solid/Tetrahedron.h"
		#include "shapes/regular_solid/TruncatedTetrahedron.h"
		#include "shapes/regular_solid/Cuboctahedron.h"
		#include "shapes/regular_solid/TruncatedCube.h"
		#include "shapes/regular_solid/TruncatedOctahedron.h"
		#include "shapes/regular_solid/Rhombicuboctahedron.h"
		#include "shapes/regular_solid/TruncatedCuboctahedron.h"
		#include "shapes/regular_solid/SnubCube.h"
		#include "shapes/regular_solid/Icosidodecahedron.h"
		#include "shapes/regular_solid/TruncatedDodecahedron.h"
		#include "shapes/regular_solid/TruncatedIcosahedron.h"
		#include "shapes/regular_solid/Rhombicosidodecahedron.h"
		#include "shapes/regular_solid/TruncatedIcosidodecahedron.h"
		#include "shapes/regular_solid/SnubDodecahedron.h"
		#include "shapes/regular_solid/CubeToTetrahedron.h"
		#include "shapes/Ellipsoid.h"
	#endif
#endif

#if RSA_SPATIAL_DIMENSION == 1
	#if RSA_ANGULAR_DIMENSION == 1
		#include "shapes/Ellipse1Dim.h"

#endif
#elif RSA_SPATIAL_DIMENSION == 2
	// 2D shapes with angular dimension
	#if RSA_ANGULAR_DIMENSION == 0
        #include "shapes/Ellipse.h"
        #include "shapes/Polydisk.h"
        #include "shapes/ordered_shape.2.1/OrderedShape2_1.h"
    #else
		#include "shapes/spherocylinder/SpheroCylinder2D.h"
		#include "shapes/Ellipse.h"
		#include "shapes/Rectangle.h"
		#include "shapes/Polydisk.h"
		#include "shapes/polygon/Polygon.h"
		#include "shapes/polygon/SBPolygon.h"
		#include "shapes/polygon/HBPolygon.h"
        #include "shapes/polygon/Triangle.h"
    #endif
#endif

struct UnknownShapeException : public std::domain_error {
    UnknownShapeException() : domain_error("") { }
};

void ShapeFactory::initShapeClass(const std::string &sClass, const std::string &attr) {
    std::ostringstream error;
    error << "[ShapeFactory::initShapeClass] ";
    try {
        ShapeFactory::initShapeClass0(sClass, attr);
        return;
    } catch (UnknownShapeException &e) {
        error << "Unknown shape: " << sClass << " or wrong dimensions: " << RSA_SPATIAL_DIMENSION << ", " << RSA_ANGULAR_DIMENSION;
    } catch (ValidationException &e) {
        error << sClass << " attributes:" << std::endl;
        error << "    " << attr << std::endl;
        error << "malformed. Reason:" << std::endl;
        error << "    " << e.what();
    }
    die(error.str());
}

void ShapeFactory::initShapeClass0(const std::string &sClass, const std::string &attr) {
    // Shapes for any dimension, without angular dimensions
    #if RSA_ANGULAR_DIMENSION == 0
        if (sClass == "Sphere") {
            Sphere<RSA_SPATIAL_DIMENSION>::initClass(attr);
            return;
        } else if (sClass == "Spherocylinder") {
             Spherocylinder<RSA_SPATIAL_DIMENSION>::initClass(attr);
             return;
        } else if (sClass == "OrientedCuboid") {
            OrientedCuboid<RSA_SPATIAL_DIMENSION>::initClass(attr);
            return;
        } else if (sClass == "OrientedPolydisk") {
            OrderedShape2_1<Polydisk>::initClass(attr);
            return;
        } else if (sClass == "OrientedEllipse") {
            OrderedShape2_1<Ellipse>::initClass(attr);
            return;
        }

		#if RSA_SPATIAL_DIMENSION == 3
        	if (sClass == "Cuboid") {
        		Cuboid::initClass(attr);
        		return;
        	} else if (sClass == "Tetrahedron") {
        		RegularSolid<Tetrahedron>::initClass(attr);
        		return;
        	} else if (sClass == "Cube") {
        		Cuboid::initClass("3 1 1 1");
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
        	} else if (sClass == "TruncatedTetrahedron") {
        		RegularSolid<TruncatedTetrahedron>::initClass(attr);
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
        	} else if (sClass == "TruncatedCuboctahedron") {
        		RegularSolid<TruncatedCuboctahedron>::initClass(attr);
        		return;
        	} else if (sClass == "SnubCube") {
        		RegularSolid<SnubCube>::initClass(attr);
        		return;
        	} else if (sClass == "Icosidodecahedron") {
        		RegularSolid<Icosidodecahedron>::initClass(attr);
        		return;
        	} else if (sClass == "TruncatedDodecahedron") {
        		RegularSolid<TruncatedDodecahedron>::initClass(attr);
        		return;
        	} else if (sClass == "TruncatedIcosahedron") {
        		RegularSolid<TruncatedIcosahedron>::initClass(attr);
        		return;
        	} else if (sClass == "Rhombicosidodecahedron") {
        		RegularSolid<Rhombicosidodecahedron>::initClass(attr);
        		return;
        	} else if (sClass == "TruncatedIcosidodecahedron") {
        		RegularSolid<TruncatedIcosidodecahedron>::initClass(attr);
        		return;
        	} else if (sClass == "SnubDodecahedron") {
        		RegularSolid<SnubDodecahedron>::initClass(attr);
        		return;
        	} else if (sClass == "Ellipsoid") {
        		Ellipsoid::initClass(attr);
        		return;
        	} else if (sClass == "CubeToTetrahedron") {
        		RegularSolid<CubeToTetrahedron>::initClass(attr);
        		return;
        	}
		#endif
	#endif
    // Shapes of specific dimensions
    #if RSA_SPATIAL_DIMENSION == 1
        // 1D shapes with angular dimension
        #if RSA_ANGULAR_DIMENSION == 1
            if (sClass == "Ellipse1Dim") {
                Ellipse1Dim::initClass(attr);
                return;
            }
        #endif

    #elif RSA_SPATIAL_DIMENSION == 2
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
            } else if (sClass == "Triangle") {
                Triangle::initClass(attr);
                return;
            } else if (sClass == "Polygon") {
                Polygon::initClass(attr);
                return;
            } else if (sClass == "Polydisk") {
              	Polydisk::initClass(attr);
              	return;
            }
        #endif

	#endif

	throw UnknownShapeException();
}

RSAShape *ShapeFactory::createShape(RND *rnd) {
    // Fetch appropriate function from Shape
    auto createShapeImpl = RSAShape::getCreateShapeImpl();
    return createShapeImpl(rnd);
}
