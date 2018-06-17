//
// Created by PKua on 26.04.18.
//

#ifndef RSA3D_OVERLAPSTRATEGYSHAPE_H
#define RSA3D_OVERLAPSTRATEGYSHAPE_H

#include "../BoundaryConditions.h"
#include "Shape.h"
#include "OverlapStrategy.h"

template <unsigned short SD, unsigned short AD>
class OverlapStrategyShape {
public:
    virtual ~OverlapStrategyShape() = default;
    virtual std::vector<std::string> getSupportedStrategies() const = 0;
    virtual OverlapStrategy<SD, AD> *createStrategy(const std::string &name) const = 0;
};

/**
 * @brief A shortcut for OverlapStrategyShape template specialization with curently used parameters.
 */
using RSAOverlapStrategyShape = OverlapStrategyShape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>;

#endif //RSA3D_OVERLAPSTRATEGYSHAPE_H
