/*
 * main.cpp
 *
 *  Created on: 20.04.2018
 *      Author: ciesla
 */

#include "../rsa3d/shape/shapes/Rectangle.h"
#include "../rsa3d/FreeBC.h"
#include  <iostream>

namespace rectanglecase
{
	int main(int argc, char **argv){

		Rectangle::initClass("1");
		FreeBC<2> bc;

		Rectangle r1, r2;
		double a1 = 3.11356;
		r1.translate(Vector<2>{{6.21535, 4.96956}});
		r1.rotate({{a1}});

		double a2 = 0.879741;
		r2.translate(Vector<2>{{7.01689, 5.90029}});
		r2.rotate({{a2}});

		std::cout << "overlap :" << r1.overlap(&bc, &r2) << ", " << r2.overlap(&bc, &r1) << std::endl;

		std::cout << "exclusion zones: " << r1.pointInside(&bc, r2.getPosition(), r2.getAngle() - 0.0001, r2.getAngle() + 0.0001) << ", " << r2.pointInside(&bc, r1.getPosition(), r1.getAngle() - 0.0001, r1.getAngle() + 0.0001) << std::endl;

		return 1;
	}
}
