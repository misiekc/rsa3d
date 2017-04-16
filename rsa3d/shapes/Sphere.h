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

	static constexpr double g20 = 1.0; // Gamma(2.0)
	static constexpr double g15 = 0.5*sqrt(M_PI); // Gamma(1.5)

	double r;

	static double gamma(int d);
	static double volume(int d);

public:
	Sphere();
	virtual ~Sphere();

	static void init(char* args);

	double getNeighbourListCellSize();
	double getVoxelSize();
	int overlap(BoundaryConditions *bc, Shape *s);
	double getVolume();
	int pointInside(BoundaryConditions *bc, double* da);
};

#endif /* SHAPES_SPHERE_H_ */
