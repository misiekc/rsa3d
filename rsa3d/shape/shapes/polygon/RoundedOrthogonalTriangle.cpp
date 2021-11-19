//
// Created by ciesla on 9.11.2021.
//

#include "RoundedOrthogonalTriangle.h"
#include "../../../utils/Assertions.h"
#include <sstream>
#include <math.h>

/**
 * OrthogonalTriangle is parameterized as follows.
 * the base of the triangle is equal to a,
 * but the total height of the rounded triangle is constant and has a constant width of 2: 2 = 2*radius + h => h = 2*(1 - radius)
 * and therefore radius < 1.0
 * the width of the triangle is w but the total width is w + 2*radius
 * The width to height ratio is thus ratio = (w+2*radius)/2 => w = 2*(ratio - radius)
 * Therefore ratio > radius
 */
void RoundedOrthogonalTriangle::initClass(const std::string &args) {
	std::istringstream in(args);

    ValidateMsg(in, "Expected 2 doubles as radius and side to base ratio");
	in >> RoundedPolygon::radius;
    ValidateMsg(RoundedPolygon::radius<1.0, "The radius has to be not greater than one");
	double ratio, w, h;
	w = 2*(1.0 - RoundedPolygon::radius);
	in >> ratio;
    ValidateMsg(ratio>RoundedPolygon::radius, "Width to height ratio should larger than radius");
    // length of straight interval
    h = 2*(ratio - RoundedPolygon::radius);

    std::ostringstream roundedPolygonAttributes;
    roundedPolygonAttributes.precision(std::numeric_limits< double >::max_digits10);

    roundedPolygonAttributes << radius << " " << 3 << " xy "
    		<< -0.5*w << " " << -0.5*h << " "
			<< -0.5*w << " " << 0.5*h << " "
			<< 0.5*w << " " << -0.5*h
			<< " 3 0 1 2";

    RoundedPolygon::initClass(roundedPolygonAttributes.str());
}
