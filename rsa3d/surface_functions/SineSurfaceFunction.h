//
// Created by Piotr Kubala on 17/10/2020.
//

#ifndef RSA3D_SINESURFACEFUNCTION_H
#define RSA3D_SINESURFACEFUNCTION_H

#include "../SurfaceFunction.h"
#include "../utils/Utils.h"
#include "../utils/Assertions.h"

class SineSurfaceFunction : public SurfaceFunction {
private:
    double A{};
    double k{};

    [[nodiscard]] inline double sinValue(double x) const { return this->A * std::sin(x); }
    inline static void normalizeTo2Pi(double &x0, double &x1);

public:
    SineSurfaceFunction(double A, double k) : A(A), k(k) { }

    [[nodiscard]] MinMax calculateValueRange(const RSAVector &voxelPosition, double voxelSpatialSize) const override;
    void fillInLastCoordinate(RSAVector &position) const override;
    [[nodiscard]] double getArea(double size) const override;
};


#endif //RSA3D_SINESURFACEFUNCTION_H
