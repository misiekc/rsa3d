//
// Created by PKua on 14.07.18.
//

#include "Rhombicosidodecahedron.h"

void Rhombicosidodecahedron::calculateStatic(const std::string &attr) {
    double g = goldRatio;
    double g2 = goldRatio*goldRatio;
    double g3 = goldRatio*goldRatio*goldRatio;
    double gt2 = 2*goldRatio;
    double gp2 = goldRatio + 2;

    RegularSolid<Rhombicosidodecahedron>::orientedVertices =
            {{{ 1,  1, g3}}, {{ 1, 1, -g3}}, {{1, -1,  g3}}, {{-1,  1,  g3}},
             {{-1, -1, g3}}, {{-1, 1, -g3}}, {{1, -1, -g3}}, {{-1, -1, -g3}},
             {{ 1,  g3, 1}}, {{ 1, g3, -1}}, {{1, -g3,  1}}, {{-1,  g3,  1}},
             {{-1, -g3, 1}}, {{-1, g3, -1}}, {{1, -g3, -1}}, {{-1, -g3, -1}},
             {{ g3,  1, 1}}, {{ g3, 1, -1}}, {{g3, -1,  1}}, {{-g3,  1,  1}},
             {{-g3, -1, 1}}, {{-g3, 1, -1}}, {{g3, -1, -1}}, {{-g3, -1, -1}},

             {{ g2,  g, gt2}}, {{ g2, g, -gt2}}, {{g2, -g,  gt2}}, {{-g2,  g,  gt2}},
             {{-g2, -g, gt2}}, {{-g2, g, -gt2}}, {{g2, -g, -gt2}}, {{-g2, -g, -gt2}},
             {{ g,  gt2, g2}}, {{ g, gt2, -g2}}, {{g, -gt2,  g2}}, {{-g,  gt2,  g2}},
             {{-g, -gt2, g2}}, {{-g, gt2, -g2}}, {{g, -gt2, -g2}}, {{-g, -gt2, -g2}},
             {{ gt2,  g2, g}}, {{ gt2, g2, -g}}, {{gt2, -g2,  g}}, {{-gt2,  g2,  g}},
             {{-gt2, -g2, g}}, {{-gt2, g2, -g}}, {{gt2, -g2, -g}}, {{-gt2, -g2, -g}},

             {{gp2, 0, g2}}, {{gp2, 0, -g2}}, {{-gp2, 0, g2}}, {{-gp2, 0, -g2}},
             {{0, g2, gp2}}, {{0, g2, -gp2}}, {{0, -g2, gp2}}, {{0, -g2, -gp2}},
             {{g2, gp2, 0}}, {{g2, -gp2, 0}}, {{-g2, gp2, 0}}, {{-g2, -gp2, 0}}};

    RegularSolid<Rhombicosidodecahedron>::orientedFaces =
            {{ 4,  3, 27, 50, 28}, {52, 32,  8, 11, 35}, {19, 43, 58, 45, 21}, {2, 26, 48, 24,  0},
             {36, 12, 10, 34, 54}, {20, 23, 47, 59, 44}, {33, 53, 37, 13,  9}, {5,  7, 31, 51, 29},
             {55, 38, 14, 15, 39}, {46, 22, 18, 42, 57}, {41, 56, 40, 16, 17}, {1, 25, 49, 30,  6},

             {27, 43, 19, 50}, { 2,  0,  3,  4}, { 3, 52, 35, 27}, {36, 54,  4, 28}, {44, 28, 50, 20}, {43, 35, 11, 58},
             { 0, 24, 32, 52}, {54, 34, 26,  2}, {44, 59, 12, 36}, {19, 21, 23, 20}, {45, 29, 51, 21}, {51, 31, 47, 23},
             {39, 15, 59, 47}, {15, 14, 10, 12}, {57, 42, 34, 10}, {42, 18, 48, 26}, {16, 40, 24, 48}, {40, 56,  8, 32},
             { 9, 13, 11,  8}, {13, 37, 45, 58}, {53,  5, 29, 37}, { 7, 55, 39, 31}, {38, 46, 57, 14}, {17, 16, 18, 22},
             {33,  9, 56, 41}, {49, 22, 46, 30}, { 6, 30, 38, 55}, {33, 25,  1, 53}, {41, 17, 49, 25}, { 1,  6,  7,  5},

             {27, 35, 43}, { 0, 52,  3}, {54,  2,  4}, {36, 28, 44}, {20, 50, 19},
             {21, 51, 23}, {31, 39, 47}, {15, 12, 59}, {14, 57, 10}, {42, 26, 34}, 
             {18, 16, 48}, {40, 32, 24}, {56,  9,  8}, {13, 58, 11}, {37, 29, 45}, 
             {49, 17, 22}, {30, 46, 38}, { 6, 55,  7}, {53,  1,  5}, {33, 41, 25}};
}

double Rhombicosidodecahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xAxis = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yAxis = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zAxis = std::abs(this->getOrientationMatrix().column(2) * axis);
    double g = goldRatio;
    double g2 = goldRatio*goldRatio;
    double g3 = goldRatio*goldRatio*goldRatio;
    double gt2 = 2*goldRatio;
    double gp2 = goldRatio + 2;

    double xRectHalfsize = std::max(std::max(xAxis + yAxis + g3*zAxis, g2*xAxis + g*yAxis + gt2*zAxis), gp2*xAxis + g2*zAxis);
    double yRectHalfsize = std::max(std::max(yAxis + zAxis + g3*xAxis, g2*yAxis + g*zAxis + gt2*xAxis), gp2*yAxis + g2*xAxis);
    double zRectHalfsize = std::max(std::max(zAxis + xAxis + g3*yAxis, g2*zAxis + g*xAxis + gt2*yAxis), gp2*zAxis + g2*yAxis);

    return std::max(std::max(xRectHalfsize, yRectHalfsize), zRectHalfsize) * normalizeFactor;
}