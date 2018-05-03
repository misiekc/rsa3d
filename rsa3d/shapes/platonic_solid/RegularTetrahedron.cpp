//
// Created by PKua on 21.04.18.
//

#include "RegularTetrahedron.h"


std::array<Vector<3>, 4> RegularTetrahedron::orientedVertices;
std::array<Vector<3>, 4> RegularTetrahedron::orientedFaceAxes;
std::array<Vector<3>, 6> RegularTetrahedron::orientedEdgeAxes;

void RegularTetrahedron::calculateStatic(const std::string &attr) {
    constexpr double edgeFactor = std::pow(3, 1./3) / 2;
    orientedVertices = {{Vector<3>{{1, 1, 1}} * edgeFactor,
                         Vector<3>{{1, -1, -1}} * edgeFactor,
                         Vector<3>{{-1, 1, -1}} * edgeFactor,
                         Vector<3>{{-1, -1, 1}} * edgeFactor}};

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
}

RegularTetrahedron::RegularTetrahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularTetrahedron>{orientation} {
    this->calculateVerticesAndAxes();
}

std::array<Vector<3>, 4> RegularTetrahedron::getFaceAxes() const { return this->faceAxes; }
std::array<Vector<3>, 6> RegularTetrahedron::getEdgeAxes() const { return this->edgeAxes; }
std::array<Vector<3>, 4> RegularTetrahedron::getVertices() const { return this->vertices; }

double RegularTetrahedron::projectionHalfsize(const Vector<3> &axis) const {
    throw std::runtime_error("unimplemented");
}

bool RegularTetrahedron::isSeparatingAxis(const Vector<3> &axis, const RegularTetrahedron &other,
                                          const Vector<3> &distance) const {
    interval thisInterval = this->getProjection(axis);
    interval otherInterval = other.getProjection(axis);

    return std::min(thisInterval.second, otherInterval.second) <= std::max(thisInterval.first, otherInterval.first);
}

RegularTetrahedron::interval RegularTetrahedron::getProjection(const Vector<3> & axis) const
{
    // Find enpoints of polyhedron projection (multiplied by unknown but const for axis factor)
    interval projInterval = {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    for (const auto &v : this->vertices) {
        double proj = v * axis;
        if (proj < projInterval.first)
            projInterval.first = proj;
        if (proj > projInterval.second)
            projInterval.second = proj;
    }
    return projInterval;
}

intersection::polyhedron RegularTetrahedron::getTriangles() const {
    return intersection::polyhedron{
            {{this->vertices[3], this->vertices[2], this->vertices[1]}},
            {{this->vertices[3], this->vertices[0], this->vertices[2]}},
            {{this->vertices[0], this->vertices[3], this->vertices[1]}},
            {{this->vertices[0], this->vertices[1], this->vertices[2]}}
    };
}
