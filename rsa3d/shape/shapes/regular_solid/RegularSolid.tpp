//
// Created by PKua on 21.04.18.
//

#include "../../../geometry/Vector.h"
#include "SATOverlapRS.h"
#include "TriTriOverlapRS.h"
#include <sstream>
#include <functional>
#include <iterator>
#include <fstream>

template<typename SpecificSolid>
const SATOverlapRS RegularSolid<SpecificSolid>::overlapStrategy{};

template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::initClass(const std::string &attr) {
    SpecificSolid::calculateStatic(attr);
    Assert(!orientedVertices.empty());
    RegularSolidBase::initClass(attr);

    Shape::setCreateShapeImpl([](RND *rnd) -> Shape* {
        return new SpecificSolid(Matrix<3, 3>::rotation(
                2*M_PI*rnd->nextValue(),
                std::asin(2*rnd->nextValue() - 1),
                2*M_PI*rnd->nextValue()));
    });
}

template<typename SpecificSolid>
bool RegularSolid<SpecificSolid>::overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const {
    SpecificSolid other = dynamic_cast<const SpecificSolid &>(*s);     // Make a copy
    this->applyBC(bc, &other);

    double distance2 = (this->getPosition() - other.getPosition()).norm2();
    if (distance2 > 4*circumsphereRadius*circumsphereRadius)
        return false;
    else if (distance2 < 4*insphereRadius*insphereRadius)
        return true;
    else
        return SpecificSolid::overlapStrategy.overlap(this, &other);
}

template<typename SpecificSolid>
Shape<3, 0> *RegularSolid<SpecificSolid>::clone() const {
    auto &thisSpecific = static_cast<const SpecificSolid&>(*this);
    return new SpecificSolid(thisSpecific);
}