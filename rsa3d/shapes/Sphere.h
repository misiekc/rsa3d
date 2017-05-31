/*
 * Sphere.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef SHAPES_SPHERE_H_
#define SHAPES_SPHERE_H_

#include <math.h>
#include <string>
#include "../RND.h"
#include "../Shape.h"

class Sphere : public Shape{

private:
	static double radius;
	static double neighbourListCellSize;
	static double voxelSize;
	static unsigned char staticDimension;

	static constexpr double g20 = 1.0; // Gamma(2.0)
	static constexpr double g15 = 0.5*sqrt(M_PI); // Gamma(1.5)

	double r;

	static double gamma(unsigned char d);
	static double volume(unsigned char d);

public:
	Sphere();
	Sphere(unsigned char dim);
	virtual ~Sphere();

	static void initClass(const std::string &args);
	static Shape * create(RND *rnd);

	double getNeighbourListCellSize();
	double getVoxelSize();
	int overlap(BoundaryConditions *bc, Shape *s);
	double getVolume();
	int pointInside(BoundaryConditions *bc, double* da);

	std::string toPovray();
	void store(std::ostream &f);
	void restore(std::istream &f);
};

#endif /* SHAPES_SPHERE_H_ */
