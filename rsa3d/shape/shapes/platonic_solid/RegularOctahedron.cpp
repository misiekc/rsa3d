//
// Created by PKua on 21.04.18.
//

#include "RegularOctahedron.h"

void RegularOctahedron::calculateStatic(const std::string &attr) {
    PlatonicSolid<RegularOctahedron>::orientedVertices = 
            {Vector<3>{{1,  0, 0}}, Vector<3>{{-1, 0, 0}}, Vector<3>{{0, 1,  0}},
             Vector<3>{{0, -1, 0}}, Vector<3>{{ 0, 0, 1}}, Vector<3>{{0, 0, -1}}};

    PlatonicSolid<RegularOctahedron>::orientedFaces =
            {{0, 2, 4}, {1, 5, 3}, {4, 2, 1}, {3, 5, 0}, {4, 1, 3}, {2, 0, 5}, {3, 0, 4}, {1, 2, 5}};
}

double RegularOctahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xHalfsize = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yHalfsize = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zHalfsize = std::abs(this->getOrientationMatrix().column(2) * axis);

    return std::max(std::max(xHalfsize, yHalfsize), zHalfsize) * normalizeFactor;
}

std::vector<double> RegularOctahedron::calculateOrder(const OrderCalculable *other) const {
    auto &otherOct = dynamic_cast<const RegularOctahedron&>(*other);

    double cos4sum = 0;
    double maxAbsCos = 0;
    auto faceAxes = this->getFaceAxes();
    auto otherFaceAxes = otherOct.getFaceAxes();
    for (const auto &a1 : faceAxes) {
        for (const auto &a2 : otherFaceAxes) {
            double absCos = std::abs(a1 * a2);
            if (absCos > maxAbsCos)
                maxAbsCos = absCos;

            double cos2 = absCos * absCos;
            cos4sum += cos2 * cos2;
        }
    }
    double nematicOrder = P2(maxAbsCos);
    double octahedralOrder = 27./128 * (5*cos4sum - 16);

    return {nematicOrder, octahedralOrder};
}