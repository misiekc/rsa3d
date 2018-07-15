//
// Created by PKua on 21.04.18.
//

#include "Tetrahedron.h"
#include "RegularSolid.h"
#include "UnoptimizedSATOverlap.h"
#include "TriTriOverlap.h"

const TriTriOverlap<Tetrahedron> Tetrahedron::overlapStrategy{};

void Tetrahedron::calculateStatic(const std::string &attr) {
    RegularSolid<Tetrahedron>::orientedVertices =
            {{{1, 1, 1}}, {{1, -1, -1}}, {{-1, 1, -1}}, {{-1, -1, 1}}};

    RegularSolid<Tetrahedron>::orientedFaces =
            {{3, 2, 1}, {3, 0, 2}, {0, 3, 1}, {0, 1, 2}};
}

double Tetrahedron::projectionHalfsize(const Vector<3> &axis) const {
    throw std::runtime_error("unimplemented");
}

bool Tetrahedron::pointInside(BoundaryConditions<3> *bc, const Vector<3> &position,
                                     const Orientation<0> &orientation, double orientationRange) const {
    if (RegularSolid::pointInside(bc, position, orientation, orientationRange))
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

OverlapStrategy<3, 0> *Tetrahedron::createStrategy(const std::string &name) const {
    if (name == "sat")
        return new UnoptimizedSATOverlap<Tetrahedron>();
    else
        return RegularSolid<Tetrahedron>::createStrategy(name);
}
