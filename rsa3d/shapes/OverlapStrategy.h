//
// Created by PKua on 26.04.18.
//

#include "../Shape.h"

#ifndef RSA3D_CUBOIDOVERLAPSTRATEGY_H
#define RSA3D_CUBOIDOVERLAPSTRATEGY_H

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
class OverlapStrategy {
public:
    virtual int overlap(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *first,
                        Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *second) const = 0;
};

#endif //RSA3D_OVERLAPSTRATEGY_H