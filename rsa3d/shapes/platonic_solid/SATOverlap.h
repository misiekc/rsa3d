//
// Created by PKua on 26.04.18.
//

#ifndef RSA3D_SATOVERLAP_H
#define RSA3D_SATOVERLAP_H


#include "../OverlapStrategy.h"

template <typename SpecificSolid>
class SATOverlap : public OverlapStrategy<3, 0> {
public:
    int overlap(Shape<3, 0> *first, Shape<3, 0> *second) const override;
};

#include "SATOverlap.tpp"

#endif //RSA3D_SATOVERLAP_H
