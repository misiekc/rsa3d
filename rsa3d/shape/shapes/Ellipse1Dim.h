/*
 * Ellipse.h
 *
 *  Created on: 10.08.2017
 *      Author: ciesla
 */

#ifndef SHAPES_ELLIPSE1DIM_H_
#define SHAPES_ELLIPSE1DIM_H_

#include "../Shape.h"
#include "../../RND.h"
#include "../../BoundaryConditions.h"
#include "../ConvexShape.h"

class Ellipse1Dim: public ConvexShape<1, 1> {
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
	bool pointInsideSpecialArea(BoundaryConditions<1> *bc, const Vector<1> &other, double angleFrom,
							   double angleTo) const;

protected:
    Matrix<2, 2> getAntiRotationMatrix() const;
    Matrix<2, 2> getRotationMatrix() const;
	void normalizeAngleRange(double *angleFrom, double *angleTo, double interval) const;
	double normalizeAngle(double angle, double interval) const;
    void setOrientation(const Orientation<1> &orientation) final;
    void setAngle(double angle);

public:
	static void initClass(const std::string &args);

	Shape<1, 1> *clone() const override;

	// Implicit copy c-tor and assignment
	Ellipse1Dim();
	~Ellipse1Dim() override = default;

    double getVolume() const override;
    bool overlap(BoundaryConditions<1> *bc, const Shape<1, 1> *s) const override;
    //bool pointInside(BoundaryConditions *bc, double* da) const override;
    bool pointInside(BoundaryConditions<1> *bc, const Vector<1> &da, double angleFrom, double angleTo) const;
    std::string toWolfram() const override;
	std::string toPovray() const override;
    std::string toString() const override;
    void store(std::ostream &f) const override;
	void restore(std::istream &f) override;
    bool pointInside(BoundaryConditions<1> *bc, const Vector<1> &position, const Orientation<1> &orientation,
                    double orientationRange) const final;
    double getAngle() const;
};


#endif /* SHAPES_ELLIPSE1DIM_H_ */
