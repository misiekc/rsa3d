//
// Created by PKua on 12.07.18.
//

#include "Rhombicuboctahedron.h"

void Rhombicuboctahedron::calculateStatic(const std::string &attr) {
    double xi = M_SQRT2 + 1;
    RegularSolid<Rhombicuboctahedron>::orientedVertices =
            {Vector<3>{{1, 1, xi}}, Vector<3>{{1, 1, -xi}}, Vector<3>{{1, -1, xi}}, Vector<3>{{-1, 1, xi}},
             Vector<3>{{-1, -1, xi}}, Vector<3>{{-1, 1, -xi}}, Vector<3>{{1, -1, -xi}}, Vector<3>{{-1, -1, -xi}},

             Vector<3>{{1, xi, 1}}, Vector<3>{{1, xi, -1}}, Vector<3>{{1, -xi, 1}}, Vector<3>{{-1, xi, 1}},
             Vector<3>{{-1, -xi, 1}}, Vector<3>{{-1, xi, -1}}, Vector<3>{{1, -xi, -1}}, Vector<3>{{-1, -xi, -1}},

             Vector<3>{{xi, 1, 1}}, Vector<3>{{xi, 1, -1}}, Vector<3>{{xi, -1, 1}}, Vector<3>{{-xi, 1, 1}},
             Vector<3>{{-xi, -1, 1}}, Vector<3>{{-xi, 1, -1}}, Vector<3>{{xi, -1, -1}}, Vector<3>{{-xi, -1, -1}}};

    /*RegularSolid<Rhombicuboctahedron>::orientedFaces =
            {{0, 4, 10, 14, 18, 22}, {20, 16, 12, 8, 4, 0}, {12, 17, 21, 1, 5, 8}, {1, 23, 19, 14, 10, 5},
             {2, 6, 9, 13, 16, 20}, {13, 9, 7, 3, 21, 17}, {7, 11, 15, 19, 23, 3}, {6, 2, 22, 18, 15, 11},

             {3, 23, 1, 21}, {15, 18, 14, 19}, {2, 20, 0, 22}, {4, 8, 5, 10}, {16, 13, 17, 12}, {11, 7, 9, 6}};*/
}

double Rhombicuboctahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xAxis = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yAxis = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zAxis = std::abs(this->getOrientationMatrix().column(2) * axis);

    double xi = M_SQRT2 + 1;
    double xRectHalfsize = xAxis + yAxis + xi*zAxis;
    double yRectHalfsize = xAxis + xi*yAxis + zAxis;
    double zRectHalfsize = xi*xAxis + yAxis + zAxis;

    return std::max(std::max(xRectHalfsize, yRectHalfsize), zRectHalfsize) * normalizeFactor;
}