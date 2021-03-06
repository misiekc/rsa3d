/*
 * Ellipse.h
 *
 *  Created on: 10.08.2017
 *      Author: ciesla
 */

#ifndef SHAPES_ELLIPSE_H_
#define SHAPES_ELLIPSE_H_

#include "../Shape.h"
#include "../../RND.h"
#include "../../BoundaryConditions.h"
#include "../AnisotropicShape2D.h"

class Ellipse: public AnisotropicShape2D {
private:
	static double longSemiAxis;
	static double shortSemiAxis;

	double a, b;
	double u[2], uT[2];

	void calculateU();
	double calculateF(double* r, double g) const;
	bool pointInsideUnrotated(const Vector<2> &p, double angleFrom, double angleTo) const;
	bool withinAngle(const Vector<2> &p, double angleFrom, double angleTo) const;
	bool withinAngleCheckCollision(const Vector<2> &p, double lowerAngle, double upperAngle) const;
	bool circleCollision(const Vector<2> &p, double tMin, double tMax) const;

public:
	Shape<2, 1> *clone() const override;

private:

	bool pointInsideSpecialArea(BoundaryConditions<2> *bc, const Vector<2> &other, double angleFrom,
							   double angleTo) const;
public:
	static void initClass(const std::string &args);

	// Implicit copy c-tor and assignment
	Ellipse();

    double getVolume(unsigned short dim) const override;
    void setAngle(double angle) override;
    bool overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const override;
    //bool pointInside(BoundaryConditions *bc, double* da) const override;
    bool pointInside(BoundaryConditions<2> *bc, const Vector<2> &da, double angleFrom, double angleTo) const override;
    std::string toWolfram() const override;
	std::string toPovray() const override;
    std::string toString() const override;
    void store(std::ostream &f) const override;
	void restore(std::istream &f) override;
};

#endif /* SHAPES_ELLIPSE_H_ */
