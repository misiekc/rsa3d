//
// Created by PKua on 21.04.18.
//

#include "RegularTetrahedron.h"


void RegularTetrahedron::calculateStatic(const std::string &attr) {
    auto &vertices = PlatonicSolid<RegularTetrahedron>::orientedVertices;
    auto &edgeAxes = PlatonicSolid<RegularTetrahedron>::orientedEdgeAxes;
    auto &faceAxes = PlatonicSolid<RegularTetrahedron>::orientedFaceAxes;
    
    constexpr double edgeFactor = std::pow(3, 1./3) / 2;
    vertices = {Vector<3>{{1, 1, 1}} * edgeFactor,
                         Vector<3>{{1, -1, -1}} * edgeFactor,
                         Vector<3>{{-1, 1, -1}} * edgeFactor,
                         Vector<3>{{-1, -1, 1}} * edgeFactor};

    edgeAxes = {vertices[0] - vertices[3],
                         vertices[0] - vertices[1],
                         vertices[0] - vertices[2],
                         vertices[2] - vertices[1],
                         vertices[1] - vertices[3],
                         vertices[3] - vertices[2]};

    faceAxes = {edgeAxes[3] ^ edgeAxes[4],
                         edgeAxes[2] ^ edgeAxes[0],
                         edgeAxes[0] ^ edgeAxes[1],
                         edgeAxes[1] ^ edgeAxes[2]};
}

RegularTetrahedron::RegularTetrahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularTetrahedron>{orientation} {
    this->calculateVerticesAndAxes();
}

double RegularTetrahedron::projectionHalfsize(const Vector<3> &axis) const {
    throw std::runtime_error("unimplemented");
}

bool RegularTetrahedron::isSeparatingAxis(const Vector<3> &axis, const RegularTetrahedron &other,
                                          const Vector<3> &distance) const {
    interval thisInterval = this->getProjection(axis);
    interval otherInterval = other.getProjection(axis);

    // TODO epsilon needed
    return std::min(thisInterval.second, otherInterval.second) < std::max(thisInterval.first, otherInterval.first);
}

RegularTetrahedron::interval RegularTetrahedron::getProjection(const Vector<3> & axis) const {
    // Find enpoints of polyhedron projection (multiplied by unknown but const for axis factor)
    interval projInterval = {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    auto vertices = this->getVertices();
    for (const auto &v : vertices) {
        double proj = v * axis;
        if (proj < projInterval.first)
            projInterval.first = proj;
        if (proj > projInterval.second)
            projInterval.second = proj;
    }
    return projInterval;
}

intersection::tri_polyh RegularTetrahedron::getTriangles() const {
    auto vertices = this->getVertices();
    return intersection::tri_polyh{
            {{vertices[3], vertices[2], vertices[1]}},
            {{vertices[3], vertices[0], vertices[2]}},
            {{vertices[0], vertices[3], vertices[1]}},
            {{vertices[0], vertices[1], vertices[2]}}
    };
}

bool RegularTetrahedron::pointInside(BoundaryConditions<3> *bc, const Vector<3> &position,
                                     const Orientation<0> &orientation, double orientationRange) const {
    if (PlatonicSolid::pointInside(bc, position, orientation, orientationRange))
        return true;

    // Additional spheres in vertices
    auto vertices = this->getVertices();
    for (auto vertex : vertices)
        if ((position - vertex).norm2() <= insphereRadius * insphereRadius)
            return true;

    // Additional spheres in midedges
    for (auto vi = vertices.begin(); vi != vertices.end(); vi++)
        for (auto vj = vi; vj != vertices.end(); vj++)
            if ((position - (*vi + *vj)/2.).norm2() <= insphereRadius * insphereRadius)
                return true;

    return false;
}
