/*
 * Sphere.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef SHAPES_SPHERE_H_
#define SHAPES_SPHERE_H_

#include <math.h>
#include "../Shape.h"

class Sphere : public Shape{

private:
	static double radius;
	static double neighbourListCellSize;
	static double voxelSize;
	static int dimension;

	static const double g20 = 1.0; // Gamma(2.0)
	static const double g15 = 0.5*sqrt(M_PI); // Gamma(1.5)

	double r;

	static double gamma(int d);
	static double volume(int d);

public:
	Sphere();
	virtual ~Sphere();

	static void init(char* args);

	double Sphere::getNeighbourListCellSize();
	double Sphere::getVoxelSize();
	int Sphere::overlap(BoundaryConditions *bc, Shape *s);
	double Sphere::getVolume();
	int Sphere::pointInside(BoundaryConditions *bc, double* da);
};

#endif /* SHAPES_SPHERE_H_ */
