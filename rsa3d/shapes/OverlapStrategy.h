//
// Created by PKua on 26.04.18.
//

#include "../Shape.h"

#ifndef RSA3D_OVERLAPSTRATEGY_H
#define RSA3D_OVERLAPSTRATEGY_H

template <unsigned short SD, unsigned short AD>
class OverlapStrategy {
public:
    virtual ~OverlapStrategy() = default;
    virtual int overlap(const Shape<SD, AD> *first, const Shape<SD, AD> *second) const = 0;
};

/**
 * @brief A shortcut for OverlapStrategy template specialization with curently used parameters.
 */
using RSAOverlapStrategy = OverlapStrategy<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>;

#endif //RSA3D_OVERLAPSTRATEGY_H