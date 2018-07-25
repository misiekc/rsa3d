//
// Created by PKua on 15.07.18.
//

#include "TruncatedTetrahedron.h"
#include "UnoptimizedSATOverlapRS.h"
#include "Tetrahedron.h"

const UnoptimizedSATOverlapRS TruncatedTetrahedron::overlapStrategy{};

void TruncatedTetrahedron::calculateStatic(const std::string &attr) {
    RegularSolid<TruncatedTetrahedron>::orientedVertices = 
            {{{ 3, 1,  1}}, {{ 1, 3,  1}}, {{ 1, 1,  3}}, {{-3, -1,  1}}, {{-1, -3,  1}}, {{-1, -1,  3}},
             {{-3, 1, -1}}, {{-1, 3, -1}}, {{-1, 1, -3}}, {{ 3, -1, -1}}, {{ 1, -3, -1}}, {{ 1, -1, -3}}};

    RegularSolid<TruncatedTetrahedron>::orientedFaces =
            {{10, 4, 3, 6, 8, 11}, {9, 11, 8, 7, 1, 0}, {4, 10, 9, 0, 2, 5}, {2, 1, 7, 6, 3, 5},

             {3, 4, 5}, {2, 0, 1}, {9, 10, 11}, {8, 6, 7}};
}

OverlapStrategy<3, 0> *TruncatedTetrahedron::createStrategy(const std::string &name) const {
    if (name == "sat")
        return new UnoptimizedSATOverlapRS();
    else
        return RegularSolid<TruncatedTetrahedron>::createStrategy(name);
}

std::vector<double> TruncatedTetrahedron::calculateOrder(const OrderCalculable *other) const {
    // Steal information from same oriented Tetraehdron
    Tetrahedron tetrahedron(this->getOrientationMatrix());
    return tetrahedron.calculateOrder(other);
}
