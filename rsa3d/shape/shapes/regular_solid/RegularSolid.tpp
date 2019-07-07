//
// Created by PKua on 21.04.18.
//

#include <sstream>
#include <functional>
#include <iterator>
#include <fstream>

#include "../../../geometry/Vector.h"
#include "SATOverlapRS.h"
#include "TriTriOverlapRS.h"
#include "../../../utils/Assertions.h"

template<typename SpecificSolid>
const SATOverlapRS RegularSolid<SpecificSolid>::overlapStrategy{};

template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::initClass(const std::string &attr) {
    SpecificSolid::calculateStatic(attr);
    Assert(!orientedVertices.empty());
    RegularSolidBase::initClass(attr);

    auto shapeInfo = getShapeStaticInfo();
    shapeInfo.setCreateShapeImpl([](RND *rnd) -> Shape* {
        return new SpecificSolid(Matrix<3, 3>::rotation(
                2*M_PI*rnd->nextValue(),
                std::asin(2*rnd->nextValue() - 1),
                2*M_PI*rnd->nextValue()));
    });
    Shape::setShapeStaticInfo(shapeInfo);
}

template<typename SpecificSolid>
bool RegularSolid<SpecificSolid>::overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const {
    switch (this->overlapEarlyRejection(bc, s)) {
        case EarlyRejectionResult::TRUE:      return true;
        case EarlyRejectionResult::FALSE:     return false;
        case EarlyRejectionResult::UNKNOWN:   break;
    }

    SpecificSolid other = dynamic_cast<const SpecificSolid &>(*s);     // Make a copy
    this->applyBC(bc, &other);
    return SpecificSolid::overlapStrategy.overlap(this, &other);
}

template<typename SpecificSolid>
Shape<3, 0> *RegularSolid<SpecificSolid>::clone() const {
    auto &thisSpecific = static_cast<const SpecificSolid&>(*this);
    return new SpecificSolid(thisSpecific);
}