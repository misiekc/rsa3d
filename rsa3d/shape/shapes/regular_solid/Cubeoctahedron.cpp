//
// Created by PKua on 12.07.18.
//

#include "Cubeoctahedron.h"

void Cubeoctahedron::calculateStatic(const std::string &attr) {
    RegularSolid<Cubeoctahedron>::orientedVertices =
            {Vector<3>{{1, 1, 0}}, Vector<3>{{-1, 1, 0}}, Vector<3>{{-1, -1, 0}}, Vector<3>{{1, -1, 0}},
             Vector<3>{{0, 1, 1}}, Vector<3>{{0, -1, 1}}, Vector<3>{{0, -1, -1}}, Vector<3>{{0, 1, -1}},
             Vector<3>{{1, 0, 1}}, Vector<3>{{-1, 0, 1}}, Vector<3>{{-1, 0, -1}}, Vector<3>{{1, 0, -1}}};

    /*RegularSolid<Cubeoctahedron>::orientedFaces =
            {{0, 2, 4}, {1, 5, 3}, {4, 2, 1}, {3, 5, 0}, {4, 1, 3}, {2, 0, 5}, {3, 0, 4}, {1, 2, 5}};*/
}

double Cubeoctahedron::projectionHalfsize(const Vector<3> &axis) const {
    return 0;
}
