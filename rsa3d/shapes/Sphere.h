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

template <unsigned short DIMENSION>
class Sphere : public Shape<DIMENSION, 0>{

private:
	static double radius;
	static double neighbourListCellSize;
	static double voxelSize;

	static constexpr double g20 = 1.0; // Gamma(2.0)
	static constexpr double g15 = 0.5*sqrt(M_PI); // Gamma(1.5)

	double r;

	static double gamma();
	static double volume();

public:
	Sphere();
	virtual ~Sphere();

	static void initClass(const std::string &args);
	static Shape<DIMENSION, 0> * create(RND *rnd);

	double getNeighbourListCellSize();
	double getVoxelSize();
	int overlap(BoundaryConditions *bc, Shape<DIMENSION, 0> *s);
	double getVolume();
	int pointInside(BoundaryConditions *bc, double* da);
	double minDistance(Shape<DIMENSION, 0> *s);

	std::string toPovray() const;
	void store(std::ostream &f) const;
	void restore(std::istream &f);
};

#include "Sphere.tpp"

#endif /* SHAPES_SPHERE_H_ */
