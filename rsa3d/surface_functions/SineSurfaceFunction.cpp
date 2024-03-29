//
// Created by pkua on 06.11.2021.
//

#include "SineSurfaceFunction.h"

SineSurfaceFunction::MinMax SineSurfaceFunction::calculateValueRange(const RSAVector &voxelPosition,
                                                                     double voxelSpatialSize) const
{
    if (this->k * voxelSpatialSize >= 2*M_PI)
        return {-this->A, this->A};

    double x0 = this->k * voxelPosition[0];
    double x1 = x0 + this->k * voxelSpatialSize;
    double min = std::min(this->sinValue(x0), this->sinValue(x1));
    double max = std::max(this->sinValue(x0), this->sinValue(x1));

    normalizeTo2Pi(x0, x1);
    // x0 > x1 means that x1 was in next 2pi period, so x0 is negative w.r.t. this period, hence
    if (x1 > M_PI/2 && x0 < M_PI/2)
        max = this->A;
    if (x1 > 3*M_PI/2 && x0 < 3*M_PI/2)
        min = -this->A;

    return {min, max};
}

void SineSurfaceFunction::normalizeTo2Pi(double &x0, double &x1) {
    while (x1 < 0) {
        x0 += 2 * M_PI;
        x1 += 2 * M_PI;
    }
    while (x1 >= 2*M_PI) {
        x0 -= 2 * M_PI;
        x1 -= 2 * M_PI;
    }
}

void SineSurfaceFunction::fillInLastCoordinate(RSAVector &position) const {
    Assert(RSA_SPATIAL_DIMENSION > 1);
    position[RSA_SPATIAL_DIMENSION - 1] = this->sinValue(this->k * position[0]);
}

double SineSurfaceFunction::getArea(double size) const {
    Assert(RSA_SPATIAL_DIMENSION > 1);
    Expects(size > 0);
    return 0;
}
