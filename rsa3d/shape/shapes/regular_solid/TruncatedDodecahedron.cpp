//
// Created by PKua on 13.07.18.
//

#include "TruncatedDodecahedron.h"
#include "Icosahedron.h"

namespace {
    const double g = RegularSolidBase::goldRatio;
    const double gi = 1/g;
    const double gp2 = g + 2;
    const double gt2 = 2*g;
    const double gp1 = g + 1;
}

void TruncatedDodecahedron::calculateStatic(const std::string &attr) {
    RegularSolid<TruncatedDodecahedron>::orientedVertices =
            {{{0, gi, gp2}}, {{0, gi, -gp2}}, {{0, -gi, gp2}}, {{0, -gi, -gp2}},
             {{gi, gp2, 0}}, {{gi, -gp2, 0}}, {{-gi, gp2, 0}}, {{-gi, -gp2, 0}},
             {{gp2, 0, gi}}, {{gp2, 0, -gi}}, {{-gp2, 0, gi}}, {{-gp2, 0, -gi}},

             {{ gi,  g, gt2}}, {{ gi, g, -gt2}}, {{gi, -g,  gt2}}, {{-gi,  g,  gt2}},
             {{-gi, -g, gt2}}, {{-gi, g, -gt2}}, {{gi, -g, -gt2}}, {{-gi, -g, -gt2}},
             {{ g,  gt2, gi}}, {{ g, gt2, -gi}}, {{g, -gt2,  gi}}, {{-g,  gt2,  gi}},
             {{-g, -gt2, gi}}, {{-g, gt2, -gi}}, {{g, -gt2, -gi}}, {{-g, -gt2, -gi}},
             {{  gt2, gi, g}}, {{ gt2, gi, -g}}, {{gt2, -gi,  g}}, {{-gt2,  gi,  g}},
             {{-gt2, -gi, g}}, {{-gt2, gi, -g}}, {{gt2, -gi, -g}}, {{-gt2, -gi, -g}},

             {{ g,  2, gp1}}, {{ g, 2, -gp1}}, {{g, -2,  gp1}}, {{-g,  2,  gp1}},
             {{-g, -2, gp1}}, {{-g, 2, -gp1}}, {{g, -2, -gp1}}, {{-g, -2, -gp1}},
             {{ 2,  gp1, g}}, {{ 2, gp1, -g}}, {{2, -gp1,  g}}, {{-2,  gp1,  g}},
             {{-2, -gp1, g}}, {{-2, gp1, -g}}, {{2, -gp1, -g}}, {{-2, -gp1, -g}},
             {{ gp1,  g, 2}}, {{ gp1, g, -2}}, {{gp1, -g,  2}}, {{-gp1,  g,  2}},
             {{-gp1, -g, 2}}, {{-gp1, g, -2}}, {{gp1, -g, -2}}, {{-gp1, -g, -2}}};

    RegularSolid<TruncatedDodecahedron>::orientedFaces =
            {{40, 16,  2,  0, 15, 39, 55, 31, 32, 56}, {31, 55, 47, 23, 25, 49, 57, 33, 11, 10}, 
             {39, 15, 12, 36, 44, 20,  4,  6, 23, 47}, {49, 25,  6,  4, 21, 45, 37, 13, 17, 41}, 
             {41, 17,  1,  3, 19, 43, 59, 35, 33, 57}, {10, 11, 35, 59, 51, 27, 24, 48, 56, 32}, 
             {27, 51, 43, 19, 18, 42, 50, 26,  5,  7}, {42, 18,  3,  1, 13, 37, 53, 29, 34, 58}, 
             {29, 53, 45, 21, 20, 44, 52, 28,  8,  9}, {28, 52, 36, 12,  0,  2, 14, 38, 54, 30},
             {38, 14, 16, 40, 48, 24,  7,  5, 22, 46}, {58, 34,  9,  8, 30, 54, 46, 22, 26, 50},
             
             {14,  2, 16}, { 0, 12, 15}, {39, 47, 55}, {31, 10, 32}, {40, 56, 48}, 
             {52, 44, 36}, {20, 21,  4}, { 6, 25, 23}, {49, 41, 57}, {33, 35, 11},
             {59, 43, 51}, {27,  7, 24}, { 5, 26, 22}, {46, 54, 38}, { 8, 28, 30},
             {34, 29,  9}, {42, 58, 50}, {19,  3, 18}, { 1, 17, 13}, {37, 45, 53}};
}

double TruncatedDodecahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xAxis = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yAxis = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zAxis = std::abs(this->getOrientationMatrix().column(2) * axis);

    double xRectHalfsize = std::max(std::max(gi*yAxis +gp2*zAxis, gi*xAxis + g*yAxis + gt2*zAxis), g*xAxis + 2*yAxis + gp1*zAxis);
    double yRectHalfsize = std::max(std::max(gi*zAxis +gp2*xAxis, gi*yAxis + g*zAxis + gt2*xAxis), g*yAxis + 2*zAxis + gp1*xAxis);
    double zRectHalfsize = std::max(std::max(gi*xAxis +gp2*yAxis, gi*zAxis + g*xAxis + gt2*yAxis), g*zAxis + 2*xAxis + gp1*yAxis);

    return std::max(std::max(xRectHalfsize, yRectHalfsize), zRectHalfsize) * normalizeFactor;
}

std::vector<double> TruncatedDodecahedron::calculateOrder(const OrderCalculable *other) const {
    // Steal information from same oriented Icosahedron
    Icosahedron icosahedron(this->getOrientationMatrix());
    return icosahedron.calculateOrder(other);
}
