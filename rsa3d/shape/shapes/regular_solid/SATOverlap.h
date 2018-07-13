//
// Created by PKua on 26.04.18.
//

#ifndef RSA3D_SATOVERLAP_H
#define RSA3D_SATOVERLAP_H


#include "../../OverlapStrategy.h"

template <typename SpecificSolid>
class SATOverlap : public OverlapStrategy<3, 0> {
private:
    bool isSeparatingAxis(const Vector<3> &axis, const SpecificSolid &first, const SpecificSolid &second,
                          const Vector<3> &distance) const {
        double distanceProj = std::abs(distance * axis);
        return distanceProj > first.projectionHalfsize(axis) + second.projectionHalfsize(axis);
    }

public:
    bool overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const override {
        auto firstSpecific = dynamic_cast<const SpecificSolid &>(*first);
        auto secondSpecific = dynamic_cast<const SpecificSolid &>(*second);
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
};


#endif //RSA3D_SATOVERLAP_H
