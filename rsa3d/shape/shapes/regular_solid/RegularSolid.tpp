//
// Created by PKua on 21.04.18.
//

#include "SATOverlapRS.h"
#include "TriTriOverlapRS.h"
#include "../../../utils/Assertions.h"

template<typename SpecificSolid>
std::shared_ptr<RegularSolidBase::ShapeData> RegularSolid<SpecificSolid>::staticShapeData;

template<typename SpecificSolid>
const SATOverlapRS RegularSolid<SpecificSolid>::overlapStrategy{};

template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::initClass(const std::string &attr) {
    // We need to initialize platonic solid corresponding to SpecificSolid's symmetry to calculate order parameters.
    // It will register that platonic solid in Shape<3, 0> class, but it is ok, we will overwrite it afterwards with
    // our SpecificSolid, however the platonic solid will remain accessible
    if (RegularSolid::isSpecificSolidBorrowingOrderParameters())
        RegularSolid<typename SpecificSolid::SymmetryPlatonicSolid>::initClass("");

    RegularSolid::staticShapeData = std::make_shared<ShapeData>(SpecificSolid::calculateStatic(attr));
    Assert(!RegularSolid::staticShapeData->orientedVertices.empty());
    RegularSolidBase::complementShapeData(*RegularSolid::staticShapeData);

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
bool RegularSolid<SpecificSolid>::isSpecificSolidBorrowingOrderParameters() {
    return !std::is_same_v<SpecificSolid, typename SpecificSolid::SymmetryPlatonicSolid>;
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

template<typename SpecificSolid>
std::vector<double> RegularSolid<SpecificSolid>::calculateOrder(const OrderCalculable *other) const {
    auto &otherRegularSolid = static_cast<const RegularSolid &>(*other);

    auto thisSymmetry = this->createSymmetryShape();
    auto otherSymmetry = otherRegularSolid.createSymmetryShape();

    return thisSymmetry->calculateOrder(otherSymmetry.get());
}

template<typename SpecificSolid>
std::unique_ptr<RegularSolidBase> RegularSolid<SpecificSolid>::createSymmetryShape() const {
    return std::make_unique<typename SpecificSolid::SymmetryPlatonicSolid>(this->getOrientationMatrix());
}
