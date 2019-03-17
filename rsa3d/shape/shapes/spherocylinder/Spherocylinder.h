/*
 * Spherocylinder.h
 *
 *  Created on: Mar 15, 2019
 *      Author: ciesla
 */

#ifndef SHAPE_SHAPES_SPHEROCYLINDER_H_
#define SHAPE_SHAPES_SPHEROCYLINDER_H_

#include "../../ConvexShape.h"
#include "../../OrderCalculable.h"
#include "../../OverlapStrategyShape.h"

template <unsigned short DIMENSION>
class Spherocylinder : public ConvexShape<DIMENSION, 0>, public OverlapStrategyShape<DIMENSION, 0>,  public OrderCalculable {

private:
	static constexpr double g20 = 1.0; // Gamma(2.0)
	static constexpr double g15 = 0.5*sqrt(M_PI); // Gamma(1.5)

	static double gamma(unsigned short dim);
	static double volume();

	static double length;
	static double radius;

	Matrix<DIMENSION, DIMENSION> orientation;

	Vector<DIMENSION> getEnd(short beginOrEnd) const;

    double distanceFrom(Vector<DIMENSION> *v) const;
    double distanceFrom(const Spherocylinder *s) const;

public:
    static void initClass(const std::string &attr);
    static double getLength();
    static double getRadius();

	explicit Spherocylinder(const Matrix<DIMENSION, DIMENSION> &orientation);

    bool overlap(BoundaryConditions<DIMENSION> *bc, const Shape<DIMENSION, 0> *s) const override;
    bool pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &position, const Orientation<0> &orientation, double orientationRange) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    Shape<DIMENSION, 0> *clone() const override;
	std::string toPovray() const;
	std::string toWolfram() const override;

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;

	std::vector<std::string> getSupportedStrategies() const override;
	OverlapStrategy<DIMENSION, 0> *createStrategy(const std::string &name) const override;

    double getAngle() const;
};

#include "Spherocylinder.tpp"

#endif /* SHAPE_SHAPES_SPHEROCYLINDER_H_ */
