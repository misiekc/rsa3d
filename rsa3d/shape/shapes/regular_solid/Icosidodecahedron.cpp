//
// Created by PKua on 13.07.18.
//

#include "Icosidodecahedron.h"
#include "Icosahedron.h"

namespace {
    const double g = RegularSolidBase::goldRatio;
    const double g2 = g*g;
    const double dg = 2*g;
}

void Icosidodecahedron::calculateStatic(const std::string &attr) {
    RegularSolid<Icosidodecahedron>::orientedVertices =
            {{{0, 0, dg}}, {{0, dg, 0}}, {{dg, 0, 0}}, {{0, 0, -dg}}, {{0, -dg, 0}}, {{-dg, 0, 0}},

             {{ 1,  g, g2}}, {{ 1, g, -g2}}, {{1, -g,  g2}}, {{-1,  g,  g2}},
             {{-1, -g, g2}}, {{-1, g, -g2}}, {{1, -g, -g2}}, {{-1, -g, -g2}},

             {{ g,  g2, 1}}, {{ g, g2, -1}}, {{g, -g2,  1}}, {{-g,  g2,  1}},
             {{-g, -g2, 1}}, {{-g, g2, -1}}, {{g, -g2, -1}}, {{-g, -g2, -1}},

             {{ g2,  1, g}}, {{ g2, 1, -g}}, {{g2, -1,  g}}, {{-g2,  1,  g}},
             {{-g2, -1, g}}, {{-g2, 1, -g}}, {{g2, -1, -g}}, {{-g2, -1, -g}}};

    RegularSolid<Icosidodecahedron>::orientedFaces =
            {{ 0, 9, 25, 26, 10}, { 6, 14,  1, 17,  9}, { 8, 24, 22,  6,  0}, { 2, 23, 15, 14, 22}, 
             {15, 7, 11, 19,  1}, {19, 27,  5, 25, 17}, {11,  3, 13, 29, 27}, {29, 21, 18, 26,  5}, 
             {18, 4, 16,  8, 10}, {21, 13, 12, 20,  4}, {20, 28,  2, 24, 16}, { 3,  7, 23, 28, 12},
             
             {22, 14,  6}, { 6,  9,  0}, { 8, 0, 10}, {24,  8, 16}, { 2, 22, 24}, 
             {28, 23,  2}, { 7, 15, 23}, {15, 1, 14}, { 1, 19, 17}, {17, 25,  9}, 
             { 5, 26, 25}, {26, 18, 10}, {21, 4, 18}, { 4, 20, 16}, {12, 28, 20}, 
             { 3, 11,  7}, {27, 19, 11}, {13, 3, 12}, {21, 29, 13}, {29,  5, 27}};
}

double Icosidodecahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xAxis = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yAxis = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zAxis = std::abs(this->getOrientationMatrix().column(2) * axis);

    double xRectHalfsize = std::max(dg*xAxis, xAxis + g*yAxis + g2*zAxis);
    double yRectHalfsize = std::max(dg*yAxis, yAxis + g*zAxis + g2*xAxis);
    double zRectHalfsize = std::max(dg*zAxis, zAxis + g*xAxis + g2*yAxis);

    return std::max(std::max(xRectHalfsize, yRectHalfsize), zRectHalfsize) * normalizeFactor;
}

std::vector<double> Icosidodecahedron::calculateOrder(const OrderCalculable *other) const {
    // Steal information from same oriented Icosahedron
    Icosahedron icosahedron(this->getOrientationMatrix());
    return icosahedron.calculateOrder(other);
}
