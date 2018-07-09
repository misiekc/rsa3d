//
// Created by PKua on 21.04.18.
//

#include "RegularTetrahedron.h"


void RegularTetrahedron::calculateStatic(const std::string &attr) {
    PlatonicSolid<RegularTetrahedron>::orientedVertices =
            {Vector<3>{{1, 1, 1}}, Vector<3>{{1, -1, -1}}, Vector<3>{{-1, 1, -1}}, Vector<3>{{-1, -1, 1}}};
    
    PlatonicSolid<RegularTetrahedron>::orientedFaces = 
            {{3, 2, 1}, {3, 0, 2}, {0, 3, 1}, {0, 1, 2}};
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
