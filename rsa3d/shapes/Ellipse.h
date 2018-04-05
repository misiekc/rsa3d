/*
 * Ellipse.h
 *
 *  Created on: 10.08.2017
 *      Author: ciesla
 */

#ifndef SHAPES_ELLIPSE_H_
#define SHAPES_ELLIPSE_H_

#include "../Shape.h"
#include "../RND.h"
#include "../BoundaryConditions.h"
#include "../AnisotropicShape2D.h"

class Ellipse: public AnisotropicShape2D {
private:
	static double longSemiAxis;
	static double shortSemiAxis;
	static double neighbourListCellSize;
	static double voxelSize;

	double a, b;
	double u[2], uT[2];

	void calculateU();
	double calculateF(double* r, double g) const;
	bool pointInsideUnrotated(const Vector<2> &p, double angleFrom, double angleTo) const;
	bool withinAngle(const Vector<2> &p, double angleFrom, double angleTo) const;
	bool withinAngleCheckCollision(const Vector<2> &p, double lowerAngle, double upperAngle) const;
	bool circleCollision(const Vector<2> &p, double tMin, double tMax) const;
	int pointInsideSpecialArea(BoundaryConditions *bc, double *other, double angleFrom, double angleTo) const;

public:
	static void initClass(const std::string &args);
	static Shape<2, 1> * create(RND *rnd);

	Ellipse();
	Ellipse(const Ellipse &other);
	Ellipse & operator = (const Ellipse & el);
	~Ellipse() override = default;

    double getNeighbourListCellSize() const override;
    double getVoxelSize() const override;
    double getVoxelAngularSize() const override;
    double getVolume() const override;
    void setAngle(double angle) override;
    int overlap(BoundaryConditions *bc, Shape<2, 1> *s) const override;
    int pointInside(BoundaryConditions *bc, double* da) const override;
    int pointInside(BoundaryConditions *bc, double* da, double angleFrom, double angleTo) const override;
    std::string toWolfram() const override;
	std::string toPovray() const override;
    std::string toString() const override;
    void store(std::ostream &f) const override;
	void restore(std::istream &f) override;
};

#endif /* SHAPES_ELLIPSE_H_ */
