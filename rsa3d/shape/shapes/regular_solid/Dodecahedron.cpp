//
// Created by PKua on 21.04.18.
//

#include "Dodecahedron.h"
#include "../../OrderParameters.h"

void Dodecahedron::calculateStatic(const std::string &attr) {
    RegularSolid<Dodecahedron>::orientedVertices =
            {{{ 1,  1,  1}}, {{-1,  1,  1}}, {{-1, -1,  1}}, {{1, -1,  1}},
             {{-1, -1, -1}}, {{-1,  1, -1}}, {{ 1,  1, -1}}, {{1, -1, -1}},

             {{gold, 1/gold, 0}}, {{-gold, 1/gold, 0}}, {{-gold, -1/gold, 0}}, {{gold, -1/gold, 0}},
             {{0, gold, 1/gold}}, {{0, -gold, 1/gold}}, {{0, -gold, -1/gold}}, {{0, gold, -1/gold}},
             {{1/gold, 0, gold}}, {{1/gold, 0, -gold}}, {{-1/gold, 0, -gold}}, {{-1/gold, 0, gold}}};
    
    RegularSolid<Dodecahedron>::orientedFaces =
            {{3, 11,  8, 0, 16}, {17,  6,  8, 11, 7}, {15, 12, 0,  8,  6}, {16,  0, 12,  1, 19},
             {3, 16, 19, 2, 13}, {11,  3, 13, 14, 7}, { 4, 10, 9,  5, 18}, { 2, 19,  1,  9, 10},
             {4, 14, 13, 2, 10}, {18, 17,  7, 14, 4}, {17, 18, 5, 15,  6}, { 9,  1, 12, 15,  5}};
}

double Dodecahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xHalfsize = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yHalfsize = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zHalfsize = std::abs(this->getOrientationMatrix().column(2) * axis);

    double cubeHalfsize = xHalfsize + yHalfsize + zHalfsize;
    double xRectHalfsize = xHalfsize * gold + yHalfsize / gold;
    double yRectHalfsize = yHalfsize * gold + zHalfsize / gold;
    double zRectHalfsize = zHalfsize * gold + xHalfsize / gold;

    return std::max(std::max(std::max(cubeHalfsize, xRectHalfsize), yRectHalfsize), zRectHalfsize) * normalizeFactor;
}

std::vector<double> Dodecahedron::calculateOrder(const OrderCalculable *other) const {
    auto &otherDod = dynamic_cast<const Dodecahedron&>(*other);
    auto thisFaceAxes = this->getFaceAxes();
    auto otherFaceAxes = otherDod.getFaceAxes();

    auto params = OrderParameters::nematicAndFull(thisFaceAxes, otherFaceAxes, OrderParameters::CosExp::_6);
    return {params.nematic, 25./192*(7*params.full - 36)};
}
