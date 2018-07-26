//
// Created by PKua on 14.07.18.
//

#include "TruncatedIcosahedron.h"

namespace {
    const double g = RegularSolidBase::goldRatio;
    const double gt3 = 3*g;
    const double gp2 = g + 2;
    const double gt2 = 2*g;
    const double g3 = g*g*g;
}

void TruncatedIcosahedron::calculateStatic(const std::string &attr) {
    RegularSolid<TruncatedIcosahedron>::orientedVertices =
            {{{0, 1, gt3}}, {{0, 1, -gt3}}, {{0, -1, gt3}}, {{0, -1, -gt3}},
             {{1, gt3, 0}}, {{1, -gt3, 0}}, {{-1, gt3, 0}}, {{-1, -gt3, 0}},
             {{gt3, 0, 1}}, {{gt3, 0, -1}}, {{-gt3, 0, 1}}, {{-gt3, 0, -1}},

             {{ 1,  gp2, gt2}}, {{ 1, gp2, -gt2}}, {{1, -gp2,  gt2}}, {{-1,  gp2,  gt2}},
             {{-1, -gp2, gt2}}, {{-1, gp2, -gt2}}, {{1, -gp2, -gt2}}, {{-1, -gp2, -gt2}},
             {{ gp2,  gt2, 1}}, {{ gp2, gt2, -1}}, {{gp2, -gt2,  1}}, {{-gp2,  gt2,  1}},
             {{-gp2, -gt2, 1}}, {{-gp2, gt2, -1}}, {{gp2, -gt2, -1}}, {{-gp2, -gt2, -1}},
             {{ gt2,  1, gp2}}, {{ gt2, 1, -gp2}}, {{gt2, -1,  gp2}}, {{-gt2,  1,  gp2}},
             {{-gt2, -1, gp2}}, {{-gt2, 1, -gp2}}, {{gt2, -1, -gp2}}, {{-gt2, -1, -gp2}},

             {{ g,  2, g3}}, {{ g, 2, -g3}}, {{g, -2,  g3}}, {{-g,  2,  g3}},
             {{-g, -2, g3}}, {{-g, 2, -g3}}, {{g, -2, -g3}}, {{-g, -2, -g3}},
             {{ 2,  g3, g}}, {{ 2, g3, -g}}, {{2, -g3,  g}}, {{-2,  g3,  g}},
             {{-2, -g3, g}}, {{-2, g3, -g}}, {{2, -g3, -g}}, {{-2, -g3, -g}},
             {{ g3,  g, 2}}, {{ g3, g, -2}}, {{g3, -g,  2}}, {{-g3, g ,  2}},
             {{-g3, -g, 2}}, {{-g3, g, -2}}, {{g3, -g, -2}}, {{-g3, -g, -2}}};

    RegularSolid<TruncatedIcosahedron>::orientedFaces =
            {{ 0, 39, 31, 32, 40,  2}, { 4,  6, 47, 15, 12, 44}, {57, 11, 10, 55, 23, 25}, {15, 47, 23, 55, 31, 39}, 
             {11, 59, 27, 24, 56, 10}, {32, 56, 24, 48, 16, 40}, { 4, 45, 13, 17, 49,  6}, {25, 49, 17, 41, 33, 57},
             {59, 35, 43, 19, 51, 27}, {41,  1,  3, 43, 35, 33}, {19, 18, 50,  5,  7, 51}, { 7,  5, 46, 14, 16, 48}, 
             {46, 22, 54, 30, 38, 14}, {38, 30, 28, 36,  0,  2}, {28, 52, 20, 44, 12, 36}, {52,  8,  9, 53, 21, 20},
             {53, 29, 37, 13, 45, 21}, {34, 42,  3,  1, 37, 29}, {26, 50, 18, 42, 34, 58}, {54, 22, 26, 58,  9,  8}, 
             
             { 9, 58, 34, 29, 53}, {18, 19, 43,  3, 42}, { 1, 41, 17, 13, 37}, {59, 11, 57, 33, 35},
             { 7, 48, 24, 27, 51}, {56, 32, 31, 55, 10}, {23, 47,  6, 49, 25}, { 0, 36, 12, 15, 39}, 
             {44, 20, 21, 45,  4}, {30, 54,  8, 52, 28}, {46,  5, 50, 26, 22}, {40, 16, 14, 38,  2}};
}

double TruncatedIcosahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xAxis = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yAxis = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zAxis = std::abs(this->getOrientationMatrix().column(2) * axis);

    double xRectHalfsize = std::max(std::max(yAxis + gt3*zAxis, xAxis + gp2*yAxis + gt2*zAxis), g*xAxis + 2*yAxis + g3*zAxis);
    double yRectHalfsize = std::max(std::max(zAxis + gt3*xAxis, yAxis + gp2*zAxis + gt2*xAxis), g*yAxis + 2*zAxis + g3*xAxis);
    double zRectHalfsize = std::max(std::max(xAxis + gt3*yAxis, zAxis + gp2*xAxis + gt2*yAxis), g*zAxis + 2*xAxis + g3*yAxis);

    return std::max(std::max(xRectHalfsize, yRectHalfsize), zRectHalfsize) * normalizeFactor;
}