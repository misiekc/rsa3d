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
#include "../../../geometry/Vector.h"

class Polygon : public Shape<2, 1> {
private:
    static void normalizeVolume(std::istringstream &in);
    static bool pointInsidePolygon(const Vector<2> &point, const std::vector<Vector<2>> &vertiecs);

    bool voxelInsideFullAngleCheck(const Vector<2> &spatialCenter, double halfSpatialSize) const;
    bool pointInsidePushedVertices(const Vector<2> &point, double pushDistance) const;
    bool pointInsidePushedEdges(const Vector<2> &point, double pushDistance) const;

    static constexpr std::size_t INSPHERE_SEARCH_DIVISIONS = 20;
    static constexpr double INSPHERE_SEARCH_FACTOR = M_SQRT2;
    static constexpr double INSPHERE_SEARCH_PRECISION = 1e-8;

protected:
	//polar coordinates of all vertices
	static std::vector<double> vertexR;
	static std::vector<double> vertexTheta;
	static std::vector<std::pair<size_t, size_t>> segments;
	static std::vector<std::pair<size_t, size_t>> helperSegments;

    static void clearOldData();
    static void parseVertices(std::istringstream &in);
    static void parseSegments(std::istringstream &in);
    static void parseHelperSegments(std::istringstream &in);

    static Vector<2> getStaticVertexPosition(std::size_t index);
    static std::vector<Vector<2>> getStaticVerticesPositions();

    static double calculateCircumscribedCircleRadius();
    static double calculateInscribedCircleRadius(const Vector<2> &origin = Vector<2>{});
	static void centerPolygon();
	static void createStarHelperSegments();

	//calculate the area of the triangle made from the origin, vertex i, and vertex j
	static double getTriangleArea(size_t i, size_t j);

	//test if line segment from point 1 to 2 intersects with line segment from point 3 to 4
	static bool lineLineIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

	//same as above, except that endpoints 3 and 4 comes from a line in a voxel, and thus carry an uncertainty
	static bool lineVoxelIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double dx, double dtheta, double l3, double l4);

    static double segmentPointDistance2(const Vector<2> &s1, const Vector<2> &s2, const Vector<2> &point);

    virtual bool overlapComplexCheck(Vector<2> &position, double angle, Vector<2> &polposition, double polangle) const;
    virtual bool voxelInsideComplexCheck(const Vector<2> &spatialCenter, double halfSpatialSize, double angularCenter,
                                 double halfAngularSize) const;

    Vector<2> getVertexPosition(std::size_t index) const;
    std::vector<Vector<2>> getVerticesPositions() const;

    void vertexToPovray(std::size_t index, std::ostream &out) const;

#ifdef CUDA_ENABLED
	static void cuInit();
	static void cuFree();
#endif


public:
	static void initClass(const std::string &args);

    static const std::vector<double> &getVertexR() { return vertexR; }
    static const std::vector<double> &getVertexTheta() { return vertexTheta; }
    static const std::vector<std::pair<size_t, size_t>> &getSegments() { return segments; }
    static const std::vector<std::pair<size_t, size_t>> &getHelperSegments() { return helperSegments; }

	Shape<2, 1> *clone() const override;
	double getVolume(unsigned short dim) const override;

	bool overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const override;

#ifdef CUDA_ENABLED
	const Shape<2,1> * overlap(BoundaryConditions<2> *bc, std::vector<const Shape<2, 1> *> *shapes) const;
#endif

	bool voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition, const Orientation<1> &voxelOrientation,
					 double spatialSize, double angularSize) const override;
	std::string toPovray() const override;
	std::string toString() const override;
	std::string toWolfram() const override;

    static bool isPolygonConvex();

    static Vector<2> calculateOptimalOrigin();
};

#endif /* SHAPES_POLYGONS_POLYGON_H_ */
