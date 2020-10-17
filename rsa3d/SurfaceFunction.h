//
// Created by Piotr Kubala on 17/10/2020.
//

#ifndef RSA3D_SURFACEFUNCTION_H
#define RSA3D_SURFACEFUNCTION_H

class SurfaceFunction {
public:
    struct MinMax {
        double min{};
        double max{};
    };

    virtual ~SurfaceFunction() = default;

    [[nodiscard]] virtual MinMax calculateValueRange(const RSAVector &voxelPosition, double voxelSpatialSize) const = 0;
    virtual void fillInLastCoordinate(RSAVector &position) const = 0;
};

#endif //RSA3D_SURFACEFUNCTION_H
