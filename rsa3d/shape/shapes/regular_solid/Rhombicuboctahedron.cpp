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

    RegularSolid<Rhombicuboctahedron>::orientedFaces =
            {{ 2,  0,  3,  4}, { 0,  8, 11,  3}, { 8,  9, 13, 11}, {  9, 1,  5, 13}, {1,  6,   7,  5}, { 6, 14, 15,  7}, 
             {14, 10, 12, 15}, {10,  2,  4, 12}, {23, 20, 19, 21}, {15, 12, 20, 23}, {4,  3,  19, 20}, {19, 11, 13, 21},
             {21,  5,  7, 23}, {16, 18, 22, 17}, { 8, 16, 17,  9}, { 0,  2, 18, 16}, {18, 10, 14, 22}, {22,  6,  1, 17}, 
             
             {7, 15, 23}, {12, 4, 20}, {3, 11, 19}, {13, 5, 21}, {6, 22, 14}, {1, 9, 17}, {8, 0, 16}, {18, 2, 10}};
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