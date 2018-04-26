//
// Created by PKua on 26.04.18.
//

#ifndef RSA3D_OVERLAPSTRATEGYSHAPE_H
#define RSA3D_OVERLAPSTRATEGYSHAPE_H

#include "../BoundaryConditions.h"
#include "../Shape.h"
#include "OverlapStrategy.h"

template <unsigned short SD, unsigned short AD>
class OverlapStrategyShape {
public:
    virtual std::vector<std::string> getSupportedStrategiesNames() const = 0;
    virtual OverlapStrategy<SD, AD> *getStrategy(const std::string &name) const = 0;

    inline int overlap(const Shape<SD, AD> *first, const Shape<SD, AD> *second,
                const OverlapStrategy<SD, AD> &strategy) const {
        return strategy.overlap(first, second);
    };
};

#endif //RSA3D_OVERLAPSTRATEGYSHAPE_H
