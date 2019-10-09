//
// Created by PKua on 21.04.18.
//

#include "Tetrahedron.h"
#include "RegularSolid.h"
#include "UnoptimizedSATOverlapRS.h"
#include "TriTriOverlapRS.h"
#include "../../OrderParameters.h"

const UnoptimizedSATOverlapRS Tetrahedron::overlapStrategy{};

RegularSolidBase::ShapeData Tetrahedron::calculateStatic(const std::string &attr) {
    ShapeData shapeData;

    shapeData.orientedVertices =
            {{{1, 1, 1}}, {{1, -1, -1}}, {{-1, 1, -1}}, {{-1, -1, 1}}};

    shapeData.orientedFaces =
            {{3, 2, 1}, {3, 0, 2}, {0, 3, 1}, {0, 1, 2}};

    // Dual tetrahedron - (-1, -1, -1) as one of vertices; CubeToTetrahedron uses it
    if (attr == "dual") {
        std::transform(shapeData.orientedVertices.begin(), shapeData.orientedVertices.end(),
                       shapeData.orientedVertices.begin(), std::negate<>());
    }

    return shapeData;
}

OverlapStrategy<3, 0> *Tetrahedron::createStrategy(const std::string &name) const {
    if (name == "sat")
        return new UnoptimizedSATOverlapRS();
    else
        return RegularSolid<Tetrahedron>::createStrategy(name);
}

std::vector<double> Tetrahedron::calculateOrder(const OrderCalculable *other) const {
    auto &otherTetr = dynamic_cast<const Tetrahedron&>(*other);
    auto thisFaceAxes = this->getFaceAxes();
    auto otherFaceAxes = otherTetr.getFaceAxes();

    return {OrderParameters::tetrahedral(thisFaceAxes, otherFaceAxes),
            OrderParameters::nematicP1(thisFaceAxes, otherFaceAxes)};
}