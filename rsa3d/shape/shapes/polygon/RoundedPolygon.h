/*
 * RoundedPolygon.h
 *
 *  Created on: 27.03.2020
 *      Author: ciesla
 */

#ifndef SHAPES_POLYGONS_ROUNDEDPOLYGON_H_
#define SHAPES_POLYGONS_ROUNDEDPOLYGON_H_

#include <vector>
#include <utility>
#include <cstddef>

#include "../../AnisotropicShape2D.h"
#include "../../../geometry/Vector.h"

#include "Polygon.h"

class RoundedPolygon : public Polygon {

private:
	static double distance2pq(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double p, double q);
	static void gradientDistance2pq(std::array<double, 2> &gradient, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double p, double q);
	static void normalizeVolume();
	static double getArea();


protected:
	static double radius;

    static double calculateCircumscribedCircleRadius();

    static double lineLineDistance2(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
    static bool lineVoxelIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double dx, double dtheta, double l3, double l4);

	bool overlapComplexCheck(Vector<2> &position, double angle, Vector<2> &polposition, double polangle) const override;
    bool voxelInsideComplexCheck(const Vector<2> &spatialCenter, double halfSpatialSize, double angularCenter, double halfAngularSize) const override;

public:
//    RoundedPolygon();

	static void initClass(const std::string &args);

	Shape<2, 1> *clone() const override;

//	bool overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const override;

//	bool voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition, const Orientation<1> &voxelOrientation,

	double getVolume(unsigned short dim) const override;

	std::string toPovray() const override;
	std::string toWolfram() const override;
};

#endif /* SHAPES_POLYGONS_ROUNDEDPOLYGON_H_ */
