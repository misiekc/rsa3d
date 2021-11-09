//
// Created by Piotr Kubala on 17/10/2020.
//

#ifndef RSA3D_FLATSURFACEFUNCTION_H
#define RSA3D_FLATSURFACEFUNCTION_H

#include "../SurfaceFunction.h"
#include "../utils/Utils.h"
#include "../utils/Assertions.h"

class FlatSurfaceFunction : public SurfaceFunction {
public:
    [[nodiscard]] MinMax calculateValueRange(const RSAVector &voxelPosition, double voxelSpatialSize) const override {
        return {0, 0};
    }

    void fillInLastCoordinate(RSAVector &position) const override {
        Assert(RSA_SPATIAL_DIMENSION > 1);
        position[RSA_SPATIAL_DIMENSION - 1] = 0;
    }

    [[nodiscard]] double getArea(double size) const override {
        Assert(RSA_SPATIAL_DIMENSION > 1);
        Expects(size > 0);
        return std::pow(size, RSA_SPATIAL_DIMENSION - 1);
    }
};


#endif //RSA3D_FLATSURFACEFUNCTION_H
