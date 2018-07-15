//
// Created by PKua on 15.07.18.
//

#include "TruncatedTetrahedron.h"
#include "UnoptimizedSATOverlap.h"

const TriTriOverlap<TruncatedTetrahedron> TruncatedTetrahedron::overlapStrategy{};

void TruncatedTetrahedron::calculateStatic(const std::string &attr) {
    RegularSolid<TruncatedTetrahedron>::orientedVertices = 
            {{{ 3, 1,  1}}, {{ 1, 3,  1}}, {{ 1, 1,  3}}, {{-3, -1,  1}}, {{-1, -3,  1}}, {{-1, -1,  3}},
             {{-3, 1, -1}}, {{-1, 3, -1}}, {{-1, 1, -3}}, {{ 3, -1, -1}}, {{ 1, -3, -1}}, {{ 1, -1, -3}}};

    RegularSolid<TruncatedTetrahedron>::orientedFaces =
            {{10, 4, 3, 6, 8, 11}, {9, 11, 8, 7, 1, 0}, {4, 10, 9, 0, 2, 5}, {2, 1, 7, 6, 3, 5},

             {3, 4, 5}, {2, 0, 1}, {9, 10, 11}, {8, 6, 7}};
}

double TruncatedTetrahedron::projectionHalfsize(const Vector<3> &axis) const {
    throw std::runtime_error("unimplemented");
}

OverlapStrategy<3, 0> *TruncatedTetrahedron::createStrategy(const std::string &name) const {
    if (name == "sat")
        return new UnoptimizedSATOverlap<TruncatedTetrahedron>();
    else
        return RegularSolid<TruncatedTetrahedron>::createStrategy(name);
}

bool TruncatedTetrahedron::overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const {
    TruncatedTetrahedron other = dynamic_cast<const TruncatedTetrahedron&>(*s);     // Make a copy
    this->applyBC(bc, &other);
    UnoptimizedSATOverlap<TruncatedTetrahedron> satOverlap;
    return satOverlap.overlap(this, &other);
}