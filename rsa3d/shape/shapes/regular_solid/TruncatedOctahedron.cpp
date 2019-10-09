//
// Created by PKua on 12.07.18.
//

#include "TruncatedOctahedron.h"

RegularSolidBase::ShapeData TruncatedOctahedron::calculateStatic(const std::string &attr) {
    ShapeData shapeData;

    shapeData.orientedVertices =
            {{{0, 1, 2}}, {{0, 1, -2}}, {{0, -1, 2}}, {{0, -1, -2}},
             {{0, 2, 1}}, {{0, 2, -1}}, {{0, -2, 1}}, {{0, -2, -1}},

             {{1, 2, 0}}, {{1, -2, 0}}, {{-1, 2, 0}}, {{-1, -2, 0}},
             {{2, 1, 0}}, {{2, -1, 0}}, {{-2, 1, 0}}, {{-2, -1, 0}},

             {{2, 0, 1}}, {{2, 0, -1}}, {{-2, 0, 1}}, {{-2, 0, -1}},
             {{1, 0, 2}}, {{1, 0, -2}}, {{-1, 0, 2}}, {{-1, 0, -2}}};

    shapeData.orientedFaces =
            {{0, 4, 10, 14, 18, 22}, {20, 16, 12, 8,  4,  0}, {12, 17, 21,  1,  5, 8}, {1, 23, 19, 14, 10,  5},
             {2, 6,  9, 13, 16, 20}, {13,  9,  7, 3, 21, 17}, { 7, 11, 15, 19, 23, 3}, {6,  2, 22, 18, 15, 11},

             {3, 23, 1, 21}, {15, 18, 14, 19}, {2, 20, 0, 22}, {4, 8, 5, 10}, {16, 13, 17, 12}, {11, 7, 9, 6}};

    return shapeData;
}

double TruncatedOctahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xAxis = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yAxis = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zAxis = std::abs(this->getOrientationMatrix().column(2) * axis);

    double xRectHalfsize = std::max(yAxis + 2*zAxis, zAxis + 2*yAxis);
    double yRectHalfsize = std::max(zAxis + 2*xAxis, xAxis + 2*zAxis);
    double zRectHalfsize = std::max(xAxis + 2*yAxis, yAxis + 2*xAxis);

    return std::max(std::max(xRectHalfsize, yRectHalfsize), zRectHalfsize) * this->shapeData->normalizeFactor;
}