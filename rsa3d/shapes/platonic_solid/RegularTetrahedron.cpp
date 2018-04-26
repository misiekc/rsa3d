//
// Created by PKua on 21.04.18.
//

#include "RegularTetrahedron.h"


std::array<Vector<3>, 4> RegularTetrahedron::orientedVertices;
std::array<Vector<3>, 4> RegularTetrahedron::orientedFaceAxes;
std::array<Vector<3>, 6> RegularTetrahedron::orientedEdgeAxes;

void RegularTetrahedron::calculateStatic(const std::string &attr) {
    orientedVertices = {{Vector<3>{{1, 1, 1}},
                 Vector<3>{{1, -1, -1}},
                 Vector<3>{{-1, 1, -1}},
                 Vector<3>{{-1, -1, 1}} }};

    orientedEdgeAxes = {{orientedVertices[0] - orientedVertices[3],
                 orientedVertices[0] - orientedVertices[1],
                 orientedVertices[0] - orientedVertices[2],
                 orientedVertices[2] - orientedVertices[1],
                 orientedVertices[1] - orientedVertices[3],
                 orientedVertices[3] - orientedVertices[2]}};

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

RegularTetrahedron::RegularTetrahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularTetrahedron>{orientation} {
    this->calculateVertices();
    this->calculateAxes();
}

std::array<Vector<3>, 4> RegularTetrahedron::getFaceAxes() const {
    return faceAxes;
}

std::array<Vector<3>, 6> RegularTetrahedron::getEdgeAxes() const {
    return edgeAxes;
}

std::array<Vector<3>, 4> RegularTetrahedron::getVertices() const {
    return this->vertices;
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
    auto thisVertices = applyOrientation(orientedVertices);
    for (const auto &v : thisVertices) {
        double proj = v * axis;
        if (proj < projInterval.first)
            projInterval.first = proj;
        if (proj > projInterval.second)
            projInterval.second = proj;
    }
    return projInterval;
}
