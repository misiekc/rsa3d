//
// Created by PKua on 12.07.18.
//

#include "TruncatedCuboctahedron.h"
#include "Octahedron.h"

namespace {
    double xi = 1 + M_SQRT2;
    double eta = 1 + 2*M_SQRT2;
}

void TruncatedCuboctahedron::calculateStatic(const std::string &attr) {
    RegularSolid<TruncatedCuboctahedron>::orientedVertices =
            {{{ 1,  xi, eta}}, {{ 1, xi, -eta}}, {{1, -xi,  eta}}, {{-1,  xi,  eta}},
             {{-1, -xi, eta}}, {{-1, xi, -eta}}, {{1, -xi, -eta}}, {{-1, -xi, -eta}},

             {{ 1,  eta, xi}}, {{ 1, eta, -xi}}, {{1, -eta,  xi}}, {{-1,  eta,  xi}},
             {{-1, -eta, xi}}, {{-1, eta, -xi}}, {{1, -eta, -xi}}, {{-1, -eta, -xi}},

             {{ xi,  eta, 1}}, {{ xi, eta, -1}}, {{xi, -eta,  1}}, {{-xi,  eta,  1}},
             {{-xi, -eta, 1}}, {{-xi, eta, -1}}, {{xi, -eta, -1}}, {{-xi, -eta, -1}},

             {{ eta,  xi, 1}}, {{ eta, xi, -1}}, {{eta, -xi,  1}}, {{-eta,  xi,  1}},
             {{-eta, -xi, 1}}, {{-eta, xi, -1}}, {{eta, -xi, -1}}, {{-eta, -xi, -1}},

             {{ eta,  1, xi}}, {{ eta, 1, -xi}}, {{eta, -1,  xi}}, {{-eta,  1,  xi}},
             {{-eta, -1, xi}}, {{-eta, 1, -xi}}, {{eta, -1, -xi}}, {{-eta, -1, -xi}},

             {{ xi,  1, eta}}, {{ xi, 1, -eta}}, {{xi, -1,  eta}}, {{-xi,  1,  eta}},
             {{-xi, -1, eta}}, {{-xi, 1, -eta}}, {{xi, -1, -eta}}, {{-xi, -1, -eta}}};

    RegularSolid<TruncatedCuboctahedron>::orientedFaces =
            {{ 2, 42, 40,  0,  3, 43, 44,  4}, {11,  8, 16, 17,  9, 13, 21, 19}, {  1, 41, 46, 6,  7, 47, 45,  5},
             {14, 22, 18, 10, 12, 20, 23, 15}, {37, 39, 31, 28, 36, 35, 27, 29}, {34, 26, 30, 38, 33, 25, 24, 32},

             {42, 34, 32, 40}, {38, 46, 41, 33}, {47, 39, 37, 45}, {36, 44, 43, 35}, {24, 25, 17, 16}, {19, 21, 29, 27},
             {28, 31, 23, 20}, {18, 22, 30, 26}, { 0,  8, 11,  3}, { 9,  1,  5, 13}, { 6, 14, 15,  7}, {10,  2,  4, 12},

             {7, 15, 23, 31, 39, 47}, {20, 12, 4, 44, 36, 28}, {43,  3, 11, 19, 27, 35}, {29, 21, 13,  5, 45, 37},
             {1,  9, 17, 25, 33, 41}, {24, 16, 8,  0, 40, 32}, {34, 42,  2, 10, 18, 26}, {38, 30, 22, 14,  6, 46}};
}

double TruncatedCuboctahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xAxis = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yAxis = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zAxis = std::abs(this->getOrientationMatrix().column(2) * axis);

    double xRectHalfsize = std::max(xAxis + xi*yAxis + eta*zAxis, xAxis + eta*yAxis + xi*zAxis);
    double yRectHalfsize = std::max(yAxis + xi*zAxis + eta*xAxis, yAxis + eta*zAxis + xi*xAxis);
    double zRectHalfsize = std::max(zAxis + xi*xAxis + eta*yAxis, zAxis + eta*xAxis + xi*yAxis);

    return std::max(std::max(xRectHalfsize, yRectHalfsize), zRectHalfsize) * normalizeFactor;
}

std::vector<double> TruncatedCuboctahedron::calculateOrder(const OrderCalculable *other) const {
    // Steal information from same oriented Octahedron
    Octahedron octahedron(this->getOrientationMatrix());
    return octahedron.calculateOrder(other);
}
