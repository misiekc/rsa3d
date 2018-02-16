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

	double calculateF(double* r, double g);


public:
	Ellipse();
	Ellipse(const Ellipse &other);
	Ellipse & operator = (const Ellipse & el);
	~Ellipse() override = default;

	static void initClass(const std::string &args);
	static Shape<2> * create(RND *rnd);

	int overlap(BoundaryConditions *bc, Shape<2> *s) override;
	double getVolume() override;
	int pointInside(BoundaryConditions *bc, double* da) override;
	int pointInside(BoundaryConditions *bc, double* da, double angleFrom, double angleTo) override;
	double getNeighbourListCellSize() override;
	double getVoxelSize() override;

    void setAngle(double angle) override;
    std::string toWolfram() const override;

	std::string toPovray() const override;
	void store(std::ostream &f) const override;
	void restore(std::istream &f) override;

    std::string toString() override;

	bool withinExclusionZoneUnrotated(const Vector<2> &p, double lowerAngle, double upperAngle) const;

    bool withinAngle(const Vector<2> &p, double lowerAngle, double upperAngle) const;

	bool withinAngleCheckCollision(const Vector<2> &p, double lowerAngle, double upperAngle) const;

	bool testCircleEllipseCollision(const Vector<2> &p, double tMin, double tMax) const;
};

#endif /* SHAPES_ELLIPSE_H_ */
