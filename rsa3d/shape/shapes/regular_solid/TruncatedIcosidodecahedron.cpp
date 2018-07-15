//
// Created by PKua on 15.07.18.
//

#include "TruncatedIcosidodecahedron.h"

namespace {
    inline double max(double a, double b, double c, double d, double e) {
        return std::max(std::max(std::max(std::max(a, b), c), d), e);
    }
}

void TruncatedIcosidodecahedron::calculateStatic(const std::string &attr) {
    double g = goldRatio;
    double gi = 1/goldRatio;
    double gp3 = goldRatio + 3;
    double git2 = 2/goldRatio;
    double gt2p1 = 2*goldRatio + 1;
    double g2 = goldRatio*goldRatio;
    double gt3m1 = 3*goldRatio - 1;
    double gt2m1 = 2*goldRatio - 1;
    double gp2 = goldRatio + 2;
    double gt2 = 2*goldRatio;

    RegularSolid<TruncatedIcosidodecahedron>::orientedVertices =
            {{{ gi,  gi, gp3}}, {{ gi, gi, -gp3}}, {{gi, -gi,  gp3}}, {{-gi,  gi,  gp3}},
             {{-gi, -gi, gp3}}, {{-gi, gi, -gp3}}, {{gi, -gi, -gp3}}, {{-gi, -gi, -gp3}},
             {{ gi,  gp3, gi}}, {{ gi, gp3, -gi}}, {{gi, -gp3,  gi}}, {{-gi,  gp3,  gi}},
             {{-gi, -gp3, gi}}, {{-gi, gp3, -gi}}, {{gi, -gp3, -gi}}, {{-gi, -gp3, -gi}},
             {{ gp3,  gi, gi}}, {{ gp3, gi, -gi}}, {{gp3, -gi,  gi}}, {{-gp3,  gi,  gi}},
             {{-gp3, -gi, gi}}, {{-gp3, gi, -gi}}, {{gp3, -gi, -gi}}, {{-gp3, -gi, -gi}},

             {{ git2,  g, gt2p1}}, {{ git2, g, -gt2p1}}, {{git2, -g,  gt2p1}}, {{-git2,  g,  gt2p1}},
             {{-git2, -g, gt2p1}}, {{-git2, g, -gt2p1}}, {{git2, -g, -gt2p1}}, {{-git2, -g, -gt2p1}},
             {{ g,  gt2p1, git2}}, {{ g, gt2p1, -git2}}, {{g, -gt2p1,  git2}}, {{-g,  gt2p1,  git2}},
             {{-g, -gt2p1, git2}}, {{-g, gt2p1, -git2}}, {{g, -gt2p1, -git2}}, {{-g, -gt2p1, -git2}},
             {{ gt2p1,  git2, g}}, {{ gt2p1, git2, -g}}, {{gt2p1, -git2,  g}}, {{-gt2p1,  git2,  g}},
             {{-gt2p1, -git2, g}}, {{-gt2p1, git2, -g}}, {{gt2p1, -git2, -g}}, {{-gt2p1, -git2, -g}},

             {{ gi,  g2, gt3m1}}, {{ gi, g2, -gt3m1}}, {{gi, -g2,  gt3m1}}, {{-gi,  g2,  gt3m1}},
             {{-gi, -g2, gt3m1}}, {{-gi, g2, -gt3m1}}, {{gi, -g2, -gt3m1}}, {{-gi, -g2, -gt3m1}},
             {{ g2,  gt3m1, gi}}, {{ g2, gt3m1, -gi}}, {{g2, -gt3m1,  gi}}, {{-g2,  gt3m1,  gi}},
             {{-g2, -gt3m1, gi}}, {{-g2, gt3m1, -gi}}, {{g2, -gt3m1, -gi}}, {{-g2, -gt3m1, -gi}},
             {{ gt3m1,   gi, g2}}, {{gt3m1, gi, -g2}}, {{gt3m1, -gi,  g2}}, {{-gt3m1,  gi,  g2}},
             {{-gt3m1, -gi, g2}}, {{-gt3m1, gi, -g2}}, {{gt3m1, -gi, -g2}}, {{-gt3m1, -gi, -g2}},

             {{ gt2m1,  2, gp2}}, {{ gt2m1, 2, -gp2}}, {{gt2m1, -2,  gp2}}, {{-gt2m1,  2,  gp2}},
             {{-gt2m1, -2, gp2}}, {{-gt2m1, 2, -gp2}}, {{gt2m1, -2, -gp2}}, {{-gt2m1, -2, -gp2}},
             {{ 2,  gp2, gt2m1}}, {{ 2, gp2, -gt2m1}}, {{2, -gp2,  gt2m1}}, {{-2,  gp2,  gt2m1}},
             {{-2, -gp2, gt2m1}}, {{-2, gp2, -gt2m1}}, {{2, -gp2, -gt2m1}}, {{-2, -gp2, -gt2m1}},
             {{ gp2,  gt2m1, 2}}, {{ gp2, gt2m1, -2}}, {{gp2, -gt2m1,  2}}, {{-gp2,  gt2m1,  2}},
             {{-gp2, -gt2m1, 2}}, {{-gp2, gt2m1, -2}}, {{gp2, -gt2m1, -2}}, {{-gp2, -gt2m1, -2}},

             {{ g,  3, gt2}}, {{ g, 3, -gt2}}, {{g, -3,  gt2}}, {{-g,  3,  gt2}},
             {{-g, -3, gt2}}, {{-g, 3, -gt2}}, {{g, -3, -gt2}}, {{-g, -3, -gt2}},
             {{ 3,  gt2, g}}, { {3, gt2, -g}}, {{3, -gt2,  g}}, {{-3,  gt2,  g}},
             {{-3, -gt2, g}}, {{-3, gt2, -g}}, {{3, -gt2, -g}}, {{-3, -gt2, -g}},
             {{ gt2,  g, 3}}, {{  gt2, g, -3}}, {{gt2, -g, 3}}, {{-gt2,  g,  3}},
             {{-gt2, -g, 3}}, {{-gt2, g, -3}}, {{gt2, -g, -3}}, {{-gt2, -g, -3}}};

    RegularSolid<TruncatedIcosidodecahedron>::orientedFaces =
            {{ 52, 100, 84,  36,  12,  10,  34,  82,  98, 50}, {26, 74, 114,  66, 64, 112,  72, 24,   0,   2},
             { 87, 103, 55,  54, 102,  86,  38,  14,  15, 39}, {44, 20,  23,  47, 95, 111,  63, 60, 108,  92},
             {110,  94, 46,  22,  18,  42,  90, 106,  58, 62}, { 7, 31,  79, 119, 71,  69, 117, 77,  29,   5},
             {  3,  27, 75, 115,  67,  68, 116,  76,  28,  4}, {78, 30,   6,   1, 25,  73, 113, 65,  70, 118},
             { 40,  16, 17,  41,  89, 105,  57,  56, 104, 88}, {80, 32,   8,  11, 35,  83,  99, 51,  48,  96},
             { 97,  49, 53, 101,  85,  37,  13,   9,  33, 81}, {59, 61, 109,  93, 45,  21,  19, 43,  91, 107},

             { 26,  2,   4, 28,  52, 50}, {100,  76, 116,  92, 108,  84}, {60,  63, 39,  15, 12,  36}, {10,  14, 38, 62, 58,  34},
             {106, 90, 114, 74,  98, 82}, { 18,  16,  40,  64,  66,  42}, { 0,  24, 48,  51, 27,   3}, {68,  67, 43, 19, 20,  44},
             {111, 95, 119, 79, 103, 87}, { 86, 102,  78, 118,  94, 110}, {46,  70, 65,  41, 17,  22}, {88, 104, 80, 96, 72, 112},
             { 99, 83, 107, 91, 115, 75}, { 31,   7,   6,  30,  54,  55}, {85, 101, 77, 117, 93, 109}, {89, 113, 73, 97, 81, 105},
             { 35, 11,  13, 37,  61, 59}, { 56,  57,  33,   9,   8,  32}, {25,   1,  5,  29, 53,  49}, {21,  45, 69, 71, 47,  23},

             {52, 28,  76, 100}, { 84, 108,  60, 36}, {12, 15,  14, 10}, { 34, 58, 106,  82}, { 2,   0,   3,   4},
             {68, 44,  92, 116}, {111,  87,  39, 63}, {38, 86, 110, 62}, { 90, 42,  66, 114}, {94, 118,  70,  46},
             {18, 22,  17,  16}, { 40,  88, 112, 64}, {72, 96,  48, 24}, { 51, 99,  75,  27}, {91,  43,  67, 115},
             {19, 21,  23,  20}, { 47,  71, 119, 95}, {31, 55, 103, 79}, { 30, 78, 102,  54}, { 5,   1,   6,   7},
             {41, 65, 113,  89}, {104,  56,  32, 80}, {35, 59, 107, 83}, {117, 69,  45,  93}, {37,  85, 109,  61},
             { 8,  9,  13,  11}, { 57, 105,  81, 33}, {73, 25,  49, 97}, { 98, 74,  26,  50}, {53,  29,  77, 101}};
}

double TruncatedIcosidodecahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xAxis = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yAxis = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zAxis = std::abs(this->getOrientationMatrix().column(2) * axis);
    double g = goldRatio;
    double gi = 1/goldRatio;
    double gp3 = goldRatio + 3;
    double git2 = 2/goldRatio;
    double gt2p1 = 2*goldRatio + 1;
    double g2 = goldRatio*goldRatio;
    double gt3m1 = 3*goldRatio - 1;
    double gt2m1 = 2*goldRatio - 1;
    double gp2 = goldRatio + 2;
    double gt2 = 2*goldRatio;

    double xRectHalfsize = max(gi*(xAxis + yAxis) + gp3*zAxis, git2*xAxis + g*yAxis + gt2p1*zAxis, gi*xAxis + g2*yAxis + gt3m1*zAxis,
                               gt2m1*xAxis + 2*yAxis + gp2*zAxis, g*xAxis + 3*yAxis +gt2*zAxis);
    double yRectHalfsize = max(gi*(yAxis + zAxis) + gp3*xAxis, git2*yAxis + g*zAxis + gt2p1*xAxis, gi*yAxis + g2*zAxis + gt3m1*xAxis,
                               gt2m1*yAxis + 2*zAxis + gp2*xAxis, g*yAxis + 3*zAxis +gt2*xAxis);
    double zRectHalfsize = max(gi*(zAxis + xAxis) + gp3*yAxis, git2*zAxis + g*xAxis + gt2p1*yAxis, gi*zAxis + g2*xAxis + gt3m1*yAxis,
                               gt2m1*zAxis + 2*xAxis + gp2*yAxis, g*zAxis + 3*xAxis +gt2*yAxis);

    return std::max(std::max(xRectHalfsize, yRectHalfsize), zRectHalfsize) * normalizeFactor;
}