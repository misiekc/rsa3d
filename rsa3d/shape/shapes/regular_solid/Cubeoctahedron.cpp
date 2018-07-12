//
// Created by PKua on 12.07.18.
//

#include "Cubeoctahedron.h"

void Cubeoctahedron::calculateStatic(const std::string &attr) {
    RegularSolid<Cubeoctahedron>::orientedVertices =
            {Vector<3>{{1, 1, 0}}, Vector<3>{{-1, 1, 0}}, Vector<3>{{-1, -1, 0}}, Vector<3>{{1, -1, 0}},
             Vector<3>{{0, 1, 1}}, Vector<3>{{0, -1, 1}}, Vector<3>{{0, -1, -1}}, Vector<3>{{0, 1, -1}},
             Vector<3>{{1, 0, 1}}, Vector<3>{{-1, 0, 1}}, Vector<3>{{-1, 0, -1}}, Vector<3>{{1, 0, -1}}};

    RegularSolid<Cubeoctahedron>::orientedFaces =
            {{5, 8, 4, 9}, {3, 11, 0, 8}, {0, 7, 1, 4}, {1, 10, 2, 9}, {2, 6, 3, 5}, {7, 11, 6, 10},

             {1, 9, 4}, {2, 5, 9}, {5, 3, 8}, {8, 0, 4}, {7, 10, 1}, {10, 6, 2}, {11, 3, 6}, {0, 11, 7}};
}

double Cubeoctahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xHalfsize = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yHalfsize = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zHalfsize = std::abs(this->getOrientationMatrix().column(2) * axis);

    return std::max(std::max(xHalfsize + yHalfsize, yHalfsize + zHalfsize), zHalfsize + xHalfsize) * normalizeFactor;
}
