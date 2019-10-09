//
// Created by PKua on 21.04.18.
//

#include "Dodecahedron.h"
#include "../../OrderParameters.h"

RegularSolidBase::ShapeData Dodecahedron::calculateStatic(const std::string &attr) {
    double g = goldRatio;

    ShapeData shapeData;

    shapeData.orientedVertices =
            {{{ 1,  1,  1}}, {{-1,  1,  1}}, {{-1, -1,  1}}, {{1, -1,  1}},
             {{-1, -1, -1}}, {{-1,  1, -1}}, {{ 1,  1, -1}}, {{1, -1, -1}},

             {{g, 1/g, 0}}, {{-g, 1/g, 0}}, {{-g, -1/g, 0}}, {{g, -1/g, 0}},
             {{0, g, 1/g}}, {{0, -g, 1/g}}, {{0, -g, -1/g}}, {{0, g, -1/g}},
             {{1/g, 0, g}}, {{1/g, 0, -g}}, {{-1/g, 0, -g}}, {{-1/g, 0, g}}};

    shapeData.orientedFaces =
            {{3, 11,  8, 0, 16}, {17,  6,  8, 11, 7}, {15, 12, 0,  8,  6}, {16,  0, 12,  1, 19},
             {3, 16, 19, 2, 13}, {11,  3, 13, 14, 7}, { 4, 10, 9,  5, 18}, { 2, 19,  1,  9, 10},
             {4, 14, 13, 2, 10}, {18, 17,  7, 14, 4}, {17, 18, 5, 15,  6}, { 9,  1, 12, 15,  5}};

    return shapeData;
}

double Dodecahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xHalfsize = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yHalfsize = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zHalfsize = std::abs(this->getOrientationMatrix().column(2) * axis);

    double cubeHalfsize = xHalfsize + yHalfsize + zHalfsize;
    double xRectHalfsize = xHalfsize * goldRatio + yHalfsize / goldRatio;
    double yRectHalfsize = yHalfsize * goldRatio + zHalfsize / goldRatio;
    double zRectHalfsize = zHalfsize * goldRatio + xHalfsize / goldRatio;

    return std::max(std::max(std::max(cubeHalfsize, xRectHalfsize), yRectHalfsize), zRectHalfsize)
           * this->shapeData->normalizeFactor;
}

std::vector<double> Dodecahedron::calculateOrder(const OrderCalculable *other) const {
    auto &otherDod = dynamic_cast<const Dodecahedron&>(*other);
    auto thisFaceAxes = this->getFaceAxes();
    auto otherFaceAxes = otherDod.getFaceAxes();
    auto thisVertexAxes = this->getVertexAxes();
    auto otherVertexAxes = otherDod.getVertexAxes();

    return {OrderParameters::dodecahedral(thisFaceAxes, otherFaceAxes),
            OrderParameters::nematicP2(thisFaceAxes, otherFaceAxes),
            OrderParameters::nematicP2(thisVertexAxes, otherVertexAxes)};
}