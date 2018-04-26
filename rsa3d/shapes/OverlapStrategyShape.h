//
// Created by PKua on 26.04.18.
//

#ifndef RSA3D_OVERLAPSTRATEGYSHAPE_H
#define RSA3D_OVERLAPSTRATEGYSHAPE_H

#include "../BoundaryConditions.h"
#include "../Shape.h"
#include "OverlapStrategy.h"

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
class OverlapStrategyShape {
public:
    virtual std::vector<std::string> getSupportedStrategiesNames() const = 0;
    virtual OverlapStrategy<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *getStrategy(const std::string &name) const = 0;

    int overlap(BoundaryConditions *bc, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *other,
                const OverlapStrategy<SPATIAL_DIMENSION, ANGULAR_DIMENSION> &strategy);
};

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int OverlapStrategyShape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::
overlap(BoundaryConditions *bc, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *other,
        const OverlapStrategy<SPATIAL_DIMENSION, ANGULAR_DIMENSION> &strategy)
{
    return 0;
}

#endif //RSA3D_OVERLAPSTRATEGYSHAPE_H
