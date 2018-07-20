//
// Created by PKua on 20.07.18.
//

#include "SATOverlapRS.h"

bool SATOverlapRS::isSeparatingAxis(const Vector<3> &axis, const RegularSolidBase &first,
                                    const RegularSolidBase &second, const Vector<3> &distance) const {
    double distanceProj = std::abs(distance * axis);
    return distanceProj > first.projectionHalfsize(axis) + second.projectionHalfsize(axis);
}

bool SATOverlapRS::overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const {
    auto &firstSpecific = dynamic_cast<const RegularSolidBase &>(*first);
    auto &secondSpecific = dynamic_cast<const RegularSolidBase &>(*second);
    Vector<3> distance = secondSpecific.getPosition() - firstSpecific.getPosition();

    // Face axes for this
    for (const auto &face : firstSpecific.getFaceAxes())
        if (this->isSeparatingAxis(face, firstSpecific, secondSpecific, distance))
            return false;
    // Face axes for other
    for (const auto &face : secondSpecific.getFaceAxes())
        if (this->isSeparatingAxis(face, firstSpecific, secondSpecific, distance))
            return false;
    // Egde axes cross products
    for (const auto &thisEdge : firstSpecific.getEdgeAxes())
        for (const auto &otherEdge : secondSpecific.getEdgeAxes())
            if (this->isSeparatingAxis(thisEdge ^ otherEdge, firstSpecific, secondSpecific, distance))
                return false;

    return true;
}
