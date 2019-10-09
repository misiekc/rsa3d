//
// Created by PKua on 21.04.18.
//

#include "Octahedron.h"
#include "../../OrderParameters.h"

RegularSolidBase::ShapeData Octahedron::calculateStatic(const std::string &attr) {
    ShapeData shapeData;

    shapeData.orientedVertices =
            {{{1, 0, 0}}, {{-1, 0, 0}}, {{0, 1, 0}}, {{0, -1, 0}}, {{0, 0, 1}}, {{0, 0, -1}}};

    shapeData.orientedFaces =
            {{0, 2, 4}, {1, 5, 3}, {4, 2, 1}, {3, 5, 0}, {4, 1, 3}, {2, 0, 5}, {3, 0, 4}, {1, 2, 5}};

    return shapeData;
}

double Octahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xHalfsize = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yHalfsize = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zHalfsize = std::abs(this->getOrientationMatrix().column(2) * axis);

    return std::max(std::max(xHalfsize, yHalfsize), zHalfsize) * this->shapeData->normalizeFactor;
}

std::vector<double> Octahedron::calculateOrder(const OrderCalculable *other) const {
    auto &otherOct = dynamic_cast<const Octahedron&>(*other);
    auto thisFaceAxes = this->getFaceAxes();
    auto otherFaceAxes = otherOct.getFaceAxes();
    auto thisVertexAxes = this->getVertexAxes();
    auto otherVertexAxes = otherOct.getVertexAxes();

    return {OrderParameters::cubatic(thisVertexAxes, otherVertexAxes),
            OrderParameters::nematicP2(thisFaceAxes, otherFaceAxes),
            OrderParameters::nematicP2(thisVertexAxes, otherVertexAxes)};
}