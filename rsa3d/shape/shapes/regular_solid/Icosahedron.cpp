//
// Created by PKua on 21.04.18.
//

#include "Icosahedron.h"
#include "../../OrderParameters.h"

void Icosahedron::calculateStatic(const std::string &attr) {
    RegularSolid<Icosahedron>::orientedVertices =
            {{{0, 1, gold}}, {{0, -1, gold}}, {{0, -1, -gold}}, {{0, 1, -gold}},
             {{gold, 0, 1}}, {{gold, 0, -1}}, {{-gold, 0, -1}}, {{-gold, 0, 1}},
             {{1, gold, 0}}, {{-1, gold, 0}}, {{-1, -gold, 0}}, {{1, -gold, 0}}};

    RegularSolid<Icosahedron>::orientedFaces =
            {{ 4,  0,  1}, { 8,  0,  4}, {9,  0,  8}, { 7, 0,  9}, { 1,  0, 7},
             { 2,  6,  3}, { 2, 10,  6}, {2, 11, 10}, { 2, 5, 11}, { 2,  3, 5},
             {11,  4,  1}, {11,  5,  4}, {5,  8,  4}, { 5, 3,  8}, { 3,  9, 8},
             { 3,  6,  9}, { 6,  7,  9}, {6, 10,  7}, {10, 1,  7}, {10, 11, 1}};
}

double Icosahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xHalfsize = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yHalfsize = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zHalfsize = std::abs(this->getOrientationMatrix().column(2) * axis);

    double xRectHalfsize = yHalfsize + zHalfsize * gold;
    double yRectHalfsize = zHalfsize + xHalfsize * gold;
    double zRectHalfsize = xHalfsize + yHalfsize * gold;

    return std::max(std::max(xRectHalfsize, yRectHalfsize), zRectHalfsize) * normalizeFactor;
}

std::vector<double> Icosahedron::calculateOrder(const OrderCalculable *other) const {
    auto &otherIcos = dynamic_cast<const Icosahedron&>(*other);
    auto thisFaceAxes = this->getFaceAxes();
    auto otherFaceAxes = otherIcos.getFaceAxes();

    auto params = OrderParameters::nematicAndFull(thisFaceAxes, otherFaceAxes, OrderParameters::CosExp::_6);
    return {params.nematic, 243./1000*(7*params.full - 100)};
}
