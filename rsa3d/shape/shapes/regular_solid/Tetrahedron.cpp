//
// Created by PKua on 21.04.18.
//

#include "Tetrahedron.h"
#include "RegularSolid.h"
#include "UnoptimizedSATOverlapRS.h"
#include "TriTriOverlapRS.h"
#include "../../OrderParameters.h"

const UnoptimizedSATOverlapRS Tetrahedron::overlapStrategy{};

void Tetrahedron::calculateStatic(const std::string &attr) {
    RegularSolid<Tetrahedron>::orientedVertices =
            {{{1, 1, 1}}, {{1, -1, -1}}, {{-1, 1, -1}}, {{-1, -1, 1}}};

    RegularSolid<Tetrahedron>::orientedFaces =
            {{3, 2, 1}, {3, 0, 2}, {0, 3, 1}, {0, 1, 2}};
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