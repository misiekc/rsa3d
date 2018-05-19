/*
 * main.cpp
 *
 *  Created on: 20.04.2018
 *      Author: ciesla
 */

#include "../rsa3d/shapes/Rectangle.h"
#include "utility/MockBC.h"
#include  <iostream>

namespace rectanglecase
{
	int main(int argc, char **argv){

		Rectangle::initClass("1");
		BoundaryConditions *bc = new MockBC();

		Rectangle r1, r2;
		double pos1[] = {6.21535, 4.96956};
		double a1 = 3.11356;
		r1.translate(pos1);
		r1.rotate({{a1}});

		double pos2[] = {7.01689, 5.90029};
		double a2 = 0.879741;
		r2.translate(pos2);
		r2.rotate({{a2}});

		std::cout << "overlap :" << r1.overlap(bc, &r2) << ", " << r2.overlap(bc, &r1) << std::endl;

		std::cout << "exclusion zones: " << r1.pointInside(bc, r2.getPosition(), r2.getAngle() - 0.0001, r2.getAngle() + 0.0001) << ", " << r2.pointInside(bc, r1.getPosition(), r1.getAngle() - 0.0001, r1.getAngle() + 0.0001) << std::endl;



		delete bc;
		return 1;
	}
}
