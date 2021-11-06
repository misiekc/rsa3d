//
// Created by ciesla on 3.11.2021.
//

#include "RoundedIsoscelesTriangle.h"
#include "../../../utils/Assertions.h"
#include <sstream>
#include <math.h>

/**
 * IsoscelesRoundedTriangle is parameterized as follows.
 * the base of the triangle is equal to a,
 * but the whole base of the rounded triangle is constant and has a constant width of 2: 2 = 2*radius + a => a = 2*(1 - radius)
 * the arm width x divided by a is the ratio - the second parameter of the triangle,
 * thus x = a*ratio = 2*(1 - radius)*ratio. x has to be greater than a/2 thus ratio has to be larger than 0.5
 * the height of the triangle is then:
 * h = sqrt(x^2-(a/2)^2) = sqrt(a*a*ratio*ratio - 0.25*a*a) = a*sqrt(ratio^2-0.25) = 2*(1 - radius)*sqrt(ratio^2-0.25)
 *
 */
void RoundedIsoscelesTriangle::initClass(const std::string &args) {
	std::istringstream in(args);

    in >> RoundedPolygon::radius;
    double ratio;
    in >> ratio;
    ValidateMsg(in, "Expected 2 doubles as radius and side to base ratio");

    ValidateMsg(RoundedPolygon::radius<1.0, "The radius has to be not greater than one");
    double a, h;
	a = 2*(1.0 - radius);
    ValidateMsg(ratio>0.5, "Sides to base ratio should not be smaller than 0.5");
    // length of straight interval
    h = a*sqrt(ratio*ratio - 0.25);

    std::ostringstream roundedPolygonAttributes;
    roundedPolygonAttributes.precision(std::numeric_limits< double >::max_digits10);

    roundedPolygonAttributes << radius << " " << 3 << " xy "
    		<< -0.5*a << " " << -0.5*h << " "
			<< 0.0 << " " << 0.5*h << " "
			<< 0.5*a << " " << -0.5*h
			<< " 3 0 1 2";

    RoundedPolygon::initClass(roundedPolygonAttributes.str());
}
