//
// Created by Piotr Kubala on 24/10/2020.
//

#ifndef RSA3D_VIRTUALSINESURFACEFUNCTION_H
#define RSA3D_VIRTUALSINESURFACEFUNCTION_H

#include "../SurfaceFunction.h"
#include "../utils/Utils.h"

class VirtualSineSurfaceFunction : public SurfaceFunction {
private:
    static constexpr double PRECISION = 1e-12;

    double A{};
    double k{};
    double r{};

    double globalMin{};
    double globalMax{};

    static void normalizeTo2Pi(double &x0, double &x1);
    [[nodiscard]] double virtualSineValue(double x) const;
    [[nodiscard]] double calculateDiskCenterX(double tangentX) const;
    [[nodiscard]] double calculateDiskCenterY(double tangentX) const;
    double findValueBisectively(double x, double xmin, double xmax) const;

public:
    VirtualSineSurfaceFunction(double A, double k, double r);

    [[nodiscard]] MinMax calculateValueRange(const RSAVector &voxelPosition, double voxelSpatialSize) const override;
    void fillInLastCoordinate(RSAVector &position) const override;

    double calculateMinValue() const;
};


#endif //RSA3D_VIRTUALSINESURFACEFUNCTION_H
