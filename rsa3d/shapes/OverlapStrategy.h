//
// Created by PKua on 26.04.18.
//

#include "../Shape.h"

#ifndef RSA3D_CUBOIDOVERLAPSTRATEGY_H
#define RSA3D_CUBOIDOVERLAPSTRATEGY_H

template <unsigned short SD, unsigned short AD>
class OverlapStrategy {
public:
    virtual int overlap(const Shape<SD, AD> *first, const Shape<SD, AD> *second) const = 0;
};

#endif //RSA3D_OVERLAPSTRATEGY_H