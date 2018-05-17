/*
 * Polygon.h
 *
 *  Created on: 16.04.2018
 *      Author: ciesla
 */

#ifndef SHAPES_POLYGONS_POLYGON_H_
#define SHAPES_POLYGONS_POLYGON_H_

#include <vector>
#include <cstddef>

#include "../AnisotropicShape2D.h"

class Polygon : public AnisotropicShape2D{

private:
	//polar coordinates of all vertices
	//assume vertex 0 is linked to vertex 1, vertex 1 is linked to vertex 2, vertex 2 is linked to vertex 3, etc.
	//assume vertex (VertexR.size()-1) is linked to vertex 0
	static std::vector<double> vertexR;
	static std::vector<double> vertexTheta;

	static double inscribedCircleRadius;

	static double getCircumscribedCircleRadius();

	static double getInscribedCircleRadius();

	//calculate the area of the triangle made from the origin, vertex i, and vertex (i+1)
	static double getTriangleArea(size_t i);

	//test if line segment from point 1 to 2 intersects with line segment from point 3 to 4
	static bool lineLineIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

	//same as above, except that endpoints 3 and 4 comes from a line in a voxel, and thus carry an uncertainty
	static bool lineVoxelIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double dx, double dtheta, double l3, double l4);

public:

	static void initClass(const std::string &args);

	Polygon();
	~Polygon() override = default;
	RSAShape *clone() const override;
	double getVolume();

	int overlap(BoundaryConditions *bc, RSAShape *s) const;
	bool voxelInside(BoundaryConditions *bc, const double *voxelPosition, const std::array<double, RSA_ANGULAR_DIMENSION> &voxelOrientation, double spatialSize, double angularSize) const override;

	std::string toPovray() const;
};

#endif /* SHAPES_POLYGONS_POLYGON_H_ */
