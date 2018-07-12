//
// Created by PKua on 26.04.18.
//

#ifndef RSA3D_SATOVERLAP_H
#define RSA3D_SATOVERLAP_H


#include "../../OverlapStrategy.h"

template <typename SpecificSolid>
class SATOverlap : public OverlapStrategy<3, 0> {
public:
    bool overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const override {
        auto firstSpecific = dynamic_cast<const SpecificSolid &>(*first);
        auto secondSpecific = dynamic_cast<const SpecificSolid &>(*second);
        Vector<3> distance = secondSpecific.getPosition() - firstSpecific.getPosition();

        // Face axes for this
        for (const auto &face : firstSpecific.getFaceAxes())
            if (firstSpecific.isSeparatingAxis(face, secondSpecific, distance))
                return false;
        // Face axes for other
        for (const auto &face : secondSpecific.getFaceAxes())
            if (firstSpecific.isSeparatingAxis(face, secondSpecific, distance))
                return false;
        // Egde axes cross products
        for (const auto &thisEdge : firstSpecific.getEdgeAxes())
            for (const auto &otherEdge : secondSpecific.getEdgeAxes())
                if (firstSpecific.isSeparatingAxis(thisEdge ^ otherEdge, secondSpecific, distance))
                    return false;

        return true;
    }
};


#endif //RSA3D_SATOVERLAP_H
