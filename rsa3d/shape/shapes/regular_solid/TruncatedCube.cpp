//
// Created by PKua on 12.07.18.
//

#include "TruncatedCube.h"
#include "Octahedron.h"

namespace {
    const double xi = M_SQRT2 - 1;
}

void TruncatedCube::calculateStatic(const std::string &attr) {
    RegularSolid<TruncatedCube>::orientedVertices =
            {{{ xi,  1, 1}}, {{ xi, 1, -1}}, {{xi, -1,  1}}, {{-xi,  1,  1}},
             {{-xi, -1, 1}}, {{-xi, 1, -1}}, {{xi, -1, -1}}, {{-xi, -1, -1}},

             {{ 1,  xi, 1}}, {{ 1, xi, -1}}, {{1, -xi,  1}}, {{-1,  xi,  1}},
             {{-1, -xi, 1}}, {{-1, xi, -1}}, {{1, -xi, -1}}, {{-1, -xi, -1}},

             {{ 1,  1, xi}}, {{ 1, 1, -xi}}, {{1, -1,  xi}}, {{-1,  1,  xi}},
             {{-1, -1, xi}}, {{-1, 1, -xi}}, {{1, -1, -xi}}, {{-1, -1, -xi}}};

    RegularSolid<TruncatedCube>::orientedFaces =
            {{ 2, 10,  8,  0,  3, 11, 12,  4}, {22, 14, 9, 17, 16, 8, 10, 18}, {17, 1, 5, 21, 19, 3,  0, 16},
             {13, 15, 23, 20, 12, 11, 19, 21}, {23,  7, 6, 22, 18, 2,  4, 20}, { 5, 1, 9, 14,  6, 7, 15, 13},
             
             {15, 7, 23}, {14, 22, 6}, {10, 2, 18}, {12, 20, 4}, {8, 16, 0}, {9, 1, 17}, {5, 13, 21}, {19, 11, 3}};
}

double TruncatedCube::projectionHalfsize(const Vector<3> &axis) const {
    double xAxis = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yAxis = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zAxis = std::abs(this->getOrientationMatrix().column(2) * axis);

    double xRectHalfsize = xi * xAxis + yAxis + zAxis;
    double yRectHalfsize = xAxis + xi * yAxis + zAxis;
    double zRectHalfsize = xAxis + yAxis + xi * zAxis;

    return std::max(std::max(xRectHalfsize, yRectHalfsize), zRectHalfsize) * normalizeFactor;
}

std::vector<double> TruncatedCube::calculateOrder(const OrderCalculable *other) const {
    // Steal information from same oriented Octahedron
    Octahedron octahedron(this->getOrientationMatrix());
    return octahedron.calculateOrder(other);
}
