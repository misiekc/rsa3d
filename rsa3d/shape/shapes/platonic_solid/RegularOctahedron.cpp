//
// Created by PKua on 21.04.18.
//

#include "RegularOctahedron.h"

void RegularOctahedron::calculateStatic(const std::string &attr) {
    PlatonicSolid<RegularOctahedron>::orientedVertices = 
            {Vector<3>{{1, 0, 0}}, Vector<3>{{-1, 0, 0}}, Vector<3>{{0, 1, 0}}, 
             Vector<3>{{0, -1, 0}}, Vector<3>{{0, 0, 1}}, Vector<3>{{0, 0, -1}}};
    PlatonicSolid<RegularOctahedron>::orientedFaces = 
            {{0, 2, 4}, {1, 5, 3}, {4, 2, 1}, {3, 5, 0}, {4, 1, 3}, {2, 0, 5}, {3, 0, 4}, {1, 2, 5}};
}

double RegularOctahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xHalfsize = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yHalfsize = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zHalfsize = std::abs(this->getOrientationMatrix().column(2) * axis);

    return std::max(std::max(xHalfsize, yHalfsize), zHalfsize) * normalizeFactor;
}