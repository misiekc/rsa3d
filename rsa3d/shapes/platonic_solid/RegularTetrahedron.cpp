//
// Created by PKua on 21.04.18.
//

#include "RegularTetrahedron.h"


std::array<Vector<3>, 4> RegularTetrahedron::orientedFaceAxes;
std::array<Vector<3>, 6> RegularTetrahedron::orientedEdgeAxes;
std::array<Vector<3>, 4> RegularTetrahedron::vertices;

void RegularTetrahedron::calculateStatic(const std::string &attr) {
    vertices = {{Vector<3>{{1, 1, 1}},
                 Vector<3>{{1, -1, -1}},
                 Vector<3>{{-1, 1, -1}},
                 Vector<3>{{-1, -1, 1}} }};

    orientedEdgeAxes = {{vertices[0] - vertices[3],
                 vertices[0] - vertices[1],
                 vertices[0] - vertices[2],
                 vertices[2] - vertices[1],
                 vertices[1] - vertices[3],
                 vertices[3] - vertices[2]}};

    orientedFaceAxes = {{orientedEdgeAxes[3] ^ orientedEdgeAxes[4],
                 orientedEdgeAxes[2] ^ orientedEdgeAxes[0],
                 orientedEdgeAxes[0] ^ orientedEdgeAxes[1],
                 orientedEdgeAxes[1] ^ orientedEdgeAxes[2]}};

    // Normalize axes
    auto normalizer = [](const Vector<3> &v) {
        return v / v.norm();
    };
    std::transform(orientedEdgeAxes.begin(), orientedEdgeAxes.end(), orientedEdgeAxes.begin(), normalizer);
    std::transform(orientedFaceAxes.begin(), orientedFaceAxes.end(), orientedFaceAxes.begin(), normalizer);
}

double RegularTetrahedron::projectionHalfsize(const Vector<3> &axis) const {
    throw std::runtime_error("unimplemented");
}

bool RegularTetrahedron::isSeparatingAxis(const Vector<3> &axis, const RegularTetrahedron &other,
                                          const Vector<3> &distance) const {
    interval thisInterval = this->getProjection(axis);
    interval otherInterval = other.getProjection(axis);

    return std::min(thisInterval.second, otherInterval.second) >= std::max(thisInterval.first, otherInterval.first);
}

RegularTetrahedron::interval RegularTetrahedron::getProjection(const Vector<3> & axis) const
{
    interval projInterval = {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};

    // Find enpoints of polyhedron projection (multiplied by unknown but const for axis factor)
    auto thisVertices = applyOrientation(vertices);
    for (const auto &v : thisVertices) {
        double proj = v * axis;
        if (proj < projInterval.first)
            projInterval.first = proj;
        if (proj > projInterval.second)
            projInterval.second = proj;
    }
    return projInterval;
}
