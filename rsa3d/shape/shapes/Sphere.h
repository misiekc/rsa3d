/*
 * Sphere.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef SHAPES_SPHERE_H_
#define SHAPES_SPHERE_H_

#include <cmath>
#include <string>
#include "../../RND.h"
#include "../ConvexShape.h"

template <unsigned short DIMENSION>
class Sphere : public ConvexShape<DIMENSION, 0>{

private:
	static double radius;

	static constexpr double g20 = 1.0; // Gamma(2.0)
	static constexpr double g15 = 0.5*sqrt(M_PI); // Gamma(1.5)

	double r;

	static double gamma();
	static double volume();

public:
	Sphere();

	static void initClass(const std::string &args);
    static double getRadius();

    bool overlap(BoundaryConditions<DIMENSION> *bc, const Shape<DIMENSION, 0> *s) const override;
	double getVolume(unsigned short dim) const override;
	bool pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &da) const override;
	bool pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &position,
					 const std::array<double ,0> &orientation, double orientationRange) const override;
	double minDistance(const Shape<DIMENSION, 0> *s) const override;

	std::string toPovray() const;
    std::string toWolfram() const override;

    void store(std::ostream &f) const;
	void restore(std::istream &f);

	Shape<DIMENSION, 0> *clone() const override;
};

#include "Sphere.tpp"

#endif /* SHAPES_SPHERE_H_ */
