//
// Created by PKua on 15.07.18.
//

#include "SnubDodecahedron.h"
#include "UnoptimizedSATOverlapRS.h"
#include "Icosahedron.h"

const TriTriOverlapRS SnubDodecahedron::overlapStrategy{};

namespace {
    const double phi = (1 + std::sqrt(5)) / 2;
    const double x = std::cbrt((phi + std::sqrt(phi-5./27))/2) + std::cbrt((phi - std::sqrt(phi-5./27))/2);

    const double c0  = phi * std::sqrt(3 - (x*x)) / 2;
    const double c1  = x * phi * std::sqrt(3 - (x*x)) / 2;
    const double c2  = phi * std::sqrt((x - 1 - (1/x)) * phi) / 2;
    const double c3  = (x*x) * phi * std::sqrt(3 - (x*x)) / 2;
    const double c4  = x * phi * std::sqrt((x - 1 - (1/x)) * phi) / 2;
    const double c5  = phi * std::sqrt(1 - x + (phi + 1) / x) / 2;
    const double c6  = phi * std::sqrt(x - phi + 1) / 2;
    const double c7  = (x*x) * phi * std::sqrt((x - 1 - (1/x)) * phi) / 2;
    const double c8  = x * phi * std::sqrt(1 - x + (phi + 1) / x) / 2;
    const double c9  = std::sqrt((x + 2) * phi + 2) / 2;
    const double c10 = x * std::sqrt(x * (phi + 1) - phi) / 2;
    const double c11 = std::sqrt((x*x) * (2 * phi + 1) - phi) / 2;
    const double c12 = phi * std::sqrt((x*x) + x) / 2;
    const double c13 = (phi*phi) * std::sqrt(x * (x + phi) + 1) / (2 * x);
    const double c14 = phi * std::sqrt(x * (x + phi) + 1) / 2;
}

void SnubDodecahedron::calculateStatic(const std::string &attr) {
    RegularSolid<SnubDodecahedron>::orientedVertices =
            {{{c2, c1, -c14}}, {{c2, -c1, c14}}, {{-c2, c1, c14}}, {{-c2, -c1, -c14}},
             {{c1, c14, -c2}}, {{c1, -c14, c2}}, {{-c1, c14, c2}}, {{-c1, -c14, -c2}},
             {{c14, c2, -c1}}, {{c14, -c2, c1}}, {{-c14, c2, c1}}, {{-c14, -c2, -c1}},

             {{c0, c8, -c12}}, {{c0, -c8, c12}}, {{-c0, c8, c12}}, {{-c0, -c8, -c12}},
             {{c8, c12, -c0}}, {{c8, -c12, c0}}, {{-c8, c12, c0}}, {{-c8, -c12, -c0}},
             {{c12, c0, -c8}}, {{c12, -c0, c8}}, {{-c12, c0, c8}}, {{-c12, -c0, -c8}},

             {{c7, c6, -c11}}, {{c7, -c6, c11}}, {{-c7, c6, c11}}, {{-c7, -c6, -c11}},
             {{c6, c11, -c7}}, {{c6, -c11, c7}}, {{-c6, c11, c7}}, {{-c6, -c11, -c7}},
             {{c11, c7, -c6}}, {{c11, -c7, c6}}, {{-c11, c7, c6}}, {{-c11, -c7, -c6}},

             {{c3, c4, c13}}, {{-c3, -c4, c13}}, {{-c3, c4, -c13}}, {{c3, -c4, -c13}},
             {{c4, c13, c3}}, {{-c4, -c13, c3}}, {{-c4, c13, -c3}}, {{c4, -c13, -c3}},
             {{c13, c3, c4}}, {{-c13, -c3, c4}}, {{-c13, c3, -c4}}, {{c13, -c3, -c4}},

             {{c9, c5, c10}}, {{-c9, -c5, c10}}, {{-c9, c5, -c10}}, {{c9, -c5, -c10}},
             {{c5, c10, c9}}, {{-c5, -c10, c9}}, {{-c5, c10, -c9}}, {{c5, -c10, -c9}},
             {{c10, c9, c5}}, {{-c10, -c9, c5}}, {{-c10, c9, -c5}}, {{c10, -c9, -c5}}};

    RegularSolid<SnubDodecahedron>::orientedFaces =
            {{46, 10, 34, 18, 58}, {45, 11, 35, 19, 57}, {22, 49, 37,  2, 26}, {30, 14, 52, 40,  6},
             {48, 36,  1, 25, 21}, {16, 56, 44,  8, 32}, {20, 51, 39,  0, 24}, {33, 17, 59, 47,  9},
             {54, 42,  4, 28, 12}, {27, 23, 50, 38,  3}, {53, 41,  5, 29, 13}, { 7, 31, 15, 55, 43},

             {34, 10, 22}, {34, 22, 26}, {34, 26, 30}, {34, 30, 18}, {18, 30,  6},
             {18,  6, 42}, {58, 18, 42}, {58, 42, 54}, {58, 54, 50}, {58, 50, 46},
             {46, 50, 23}, {46, 23, 11}, {46, 11, 10}, {10, 11, 45}, {10, 45, 22},
             {11, 23, 35}, {35, 23, 27}, {35, 27, 31}, {35, 31, 19}, {19, 31,  7},
             {19,  7, 41}, {19, 41, 57}, {57, 41, 53}, {57, 53, 49}, {57, 49, 45},
             {45, 49, 22}, {49, 53, 37}, {37, 53, 13}, {37, 13,  1}, {37,  1,  2},
             { 2,  1, 36}, { 2, 36, 14}, { 2, 14, 26}, {26, 14, 30}, {14, 36, 52},
             {52, 36, 48}, {52, 48, 56}, {52, 56, 40}, {40, 56, 16}, {40, 16,  4},
             {40,  4,  6}, { 6,  4, 42}, { 1, 13, 25}, {25, 13, 29}, {25, 29, 33},
             {25, 33, 21}, {21, 33,  9}, {21,  9, 44}, {21, 44, 48}, {48, 44, 56},
             {44,  9,  8}, { 8,  9, 47}, { 8, 47, 20}, { 8, 20, 32}, {32, 20, 24},
             {32, 24, 28}, {32, 28, 16}, {16, 28,  4}, {28, 24, 12}, {12, 24,  0},
             {12,  0, 38}, {12, 38, 54}, {54, 38, 50}, {38,  0,  3}, { 3,  0, 39},
             { 3, 39, 15}, { 3, 15, 27}, {27, 15, 31}, {15, 39, 55}, {55, 39, 51},
             {55, 51, 59}, {55, 59, 43}, {43, 59, 17}, {43, 17,  5}, {43,  5,  7},
             { 7,  5, 41}, { 5, 17, 29}, {29, 33, 17}, {59, 51, 47}, {47, 51, 20}};
}

OverlapStrategy<3, 0> *SnubDodecahedron::createStrategy(const std::string &name) const {
    if (name == "sat")
        return new UnoptimizedSATOverlapRS();
    else
        return RegularSolid<SnubDodecahedron>::createStrategy(name);
}

std::vector<double> SnubDodecahedron::calculateOrder(const OrderCalculable *other) const {
    // Steal information from same oriented Icosahedron
    Icosahedron icosahedron(this->getOrientationMatrix());
    return icosahedron.calculateOrder(other);
}
