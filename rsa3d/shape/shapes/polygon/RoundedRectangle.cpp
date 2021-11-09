//
// Created by ciesla on 3.11.2021.
//

#include "RoundedRectangle.h"
#include "../../../utils/Assertions.h"
#include <sstream>


/**
 * Rectangle is parametrized as follows:
 * Its total height is constant and equal to 2. Thus 2 = 2*radius + h => h = 2(1-radius)
 * Its total width is equal to: 2*radius + w
 * The total width to the total height ratio is: ratio = (2*radius + w)/2 => w = 2*(ratio-radius)
 * @param args
 */
void RoundedRectangle::initClass(const std::string &args) {
	std::istringstream in(args);

    ValidateMsg(in, "Expected 2 doubles as radius and width to height ratio");
	in >> RoundedPolygon::radius;
    ValidateMsg(RoundedPolygon::radius<1.0, "The radius has to be not greater than one");
	// width to height ratio
	double ratio, h, w;
	h = 2*(1.0 - RoundedPolygon::radius);
	in >> ratio;
    ValidateMsg(ratio>=1.0, "Width to height ratio should not be smaller than 1.0");
    // length of straight interval: ratio = (w+2*radius)/(1+2*radius)
    w = 2*(ratio-RoundedPolygon::radius);

    std::ostringstream roundedPolygonAttributes;
    roundedPolygonAttributes.precision(std::numeric_limits< double >::max_digits10);

    roundedPolygonAttributes << radius << " " << 4 << " xy "
    		<< 0.5*h << " " << -0.5*w << " "
    		<< 0.5*h << " " << 0.5*w << " "
			<< -0.5*h << " " << 0.5*w << " "
			<< -0.5*h << " " << -0.5*w << " "
			<< " 4 0 1 2 3";
    RoundedPolygon::initClass(roundedPolygonAttributes.str());

    // Make angular voxel size smaller - regular polygon has an 2-fold rotational symmetry
    ShapeStaticInfo<2, 1> shapeInfo = Shape::getShapeStaticInfo();
    shapeInfo.setAngularVoxelSize(M_PI);
    Shape::setShapeStaticInfo(shapeInfo);
}
