/*
 * Polygon.h
 *
 *  Created on: 16.04.2018
 *      Author: ciesla
 */

#ifndef SHAPES_POLYGONS_POLYGON_H_
#define SHAPES_POLYGONS_POLYGON_H_

#include <vector>
#include <utility>
#include <cstddef>

#include "../../AnisotropicShape2D.h"
#include "../../../Vector.h"

class Polygon : public Shape<2, 1> {

private:
	//test if line segment from point 1 to 2 intersects with line segment from point 3 to 4
	static bool lineLineIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

	//same as above, except that endpoints 3 and 4 comes from a line in a voxel, and thus carry an uncertainty
	static bool lineVoxelIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double dx, double dtheta, double l3, double l4);

	Vector<2> getVertexPosition(std::size_t index) const;
	void vertexToPovray(std::size_t index, std::ostream &out) const;

    bool voxelInsideEasyCheck(const Vector<2> &spatialCenter, double halfSpatialSize) const;
    bool voxelInsideComplexCheck(const Vector<2> &spatialCenter, double halfSpatialSize, double angularCenter,
                                 double halfAngularSize) const;
    bool voxelInsideFullAngleCheck(const Vector<2> &spatialCenter, double halfSpatialSize) const;

    bool pointInsidePushed(const Vector<2> &point, double pushDistance) const;
    bool pointInsidePolygon(const Vector<2> &point) const;
    bool pointInsidePushedBoundary(const Vector<2> &point, double pushDistance) const;

protected:
	//polar coordinates of all vertices
	//assume vertex 0 is linked to vertex 1, vertex 1 is linked to vertex 2, vertex 2 is linked to vertex 3, etc.
	//assume vertex (VertexR.size()-1) is linked to vertex 0
	static std::vector<double> vertexR;
	static std::vector<double> vertexTheta;
	static std::vector<std::pair<size_t, size_t>> segments;
	static std::vector<std::pair<size_t, size_t>> helperSegments;

	static double inscribedCircleRadius;

	static double getCircumscribedCircleRadius();
	static double getInscribedCircleRadius();

	static void centerPolygon();
	static void createStarHelperSegments();

	//calculate the area of the triangle made from the origin, vertex i, and vertex j
	static double getTriangleArea(size_t i, size_t j);

	static float *d_vertices;
	static float *d_angles;
	static std::pair<unsigned int, unsigned int> *d_segments;
	static std::pair<unsigned int, unsigned int> *d_helperSegments;

#ifdef CUDA_ENABLED
	static void cuInit();
	static void cuFree();
#endif


public:

	static void initClass(const std::string &args);

	~Polygon() override = default;

	Shape<2, 1> *clone() const override;
	double getVolume() const override;

	bool overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const override;

#ifdef CUDA_ENABLED
	const Shape<2,1> * overlap(BoundaryConditions<2> *bc, std::vector<const Shape<2, 1> *> *shapes) const;
#endif

	bool voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition, const Orientation<1> &voxelOrientation,
					 double spatialSize, double angularSize) const override;
	std::string toPovray() const override;
	std::string toString() const override;
	std::string toWolfram() const override;
};

#endif /* SHAPES_POLYGONS_POLYGON_H_ */
