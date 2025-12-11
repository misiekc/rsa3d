/*
 * Sphere.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef SHAPES_SPHERE_H_
#define SHAPES_SPHERE_H_

#include <string>
#include "../../RND.h"
#include "../ConvexShape.h"

template <unsigned short DIMENSION>
class Sphere : public ConvexShape<DIMENSION, 0>{

private:
	static double radius;
	static double g20;
	static double g15;

	double r;

	static double gamma(unsigned short dim);
	static double volume(unsigned short dim);

public:
	Sphere();

	static void initClass(const std::string &args);
    static double getRadius();

    bool overlap(BoundaryConditions<DIMENSION> *bc, const Shape<DIMENSION, 0> *s) const override;
	double getVolume(unsigned short dim) const override;
	bool pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &da) const override;
	bool pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &position,
					 const Orientation<0> &orientation, const Orientation<0> &orientationRange) const override;
	double minDistance(const Shape<DIMENSION, 0> *s) const override;

	std::string toPovray() const override;
    std::string toWolfram() const override;

    void store(std::ostream &f) const override;
	void restore(std::istream &f) override;

	Shape<DIMENSION, 0> *clone() const override;
};

#include "Sphere.tpp"

#endif /* SHAPES_SPHERE_H_ */
