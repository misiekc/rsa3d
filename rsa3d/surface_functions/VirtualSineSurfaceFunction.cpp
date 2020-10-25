//
// Created by Piotr Kubala on 24/10/2020.
//

#include "VirtualSineSurfaceFunction.h"
#include "../utils/Assertions.h"

SurfaceFunction::MinMax VirtualSineSurfaceFunction::calculateValueRange(const RSAVector &voxelPosition,
                                                                        double voxelSpatialSize) const
{
    if (this->k * voxelSpatialSize >= 2*M_PI)
        return {this->calculateMinValue(), this->A + this->r};

    double x0 = voxelPosition[0];
    double x1 = x0 + voxelSpatialSize;
    double val0 = this->virtualSineValue(x0);
    double val1 = this->virtualSineValue(x1);

    double min = std::min(val0, val1);
    double max = std::max(val0, val1);

    x0 *= this->k;
    x1 *= this->k;
    normalizeTo2Pi(x0, x1);

    if ((x1 > M_PI/2 && x0 < M_PI/2) || (x0 < -3*M_PI/2))
        max = this->A + this->r;
    if ((x1 > 3*M_PI/2 && x0 < 3*M_PI/2) || (x0 < -M_PI/2))
        min = this->calculateMinValue();

    return {min, max};
}

void VirtualSineSurfaceFunction::fillInLastCoordinate(RSAVector &position) const {
    position[RSA_SPATIAL_DIMENSION - 1] = this->virtualSineValue(position[0]);
}

void VirtualSineSurfaceFunction::normalizeTo2Pi(double &x0, double &x1) {
    while (x1 < 0) {
        x0 += 2 * M_PI;
        x1 += 2 * M_PI;
    }
    while (x1 >= 2*M_PI) {
        x0 -= 2 * M_PI;
        x1 -= 2 * M_PI;
    }
}

double VirtualSineSurfaceFunction::virtualSineValue(double x) const {
    x *= this->k;
    while (x < 0)
        x += 2 * M_PI;
    while (x >= 2 * M_PI)
        x -= 2 * M_PI;

    double xmin{}, xmax{};
    if (x < M_PI/2) {
        xmin = -M_PI/2;
        xmax = M_PI/2;
    } else if (x < 3*M_PI/2) {
        xmin = M_PI/2;
        xmax = 3*M_PI/2;
    } else {
        xmin = 3*M_PI/2;
        xmax = 5*M_PI/2;
    }

    x /= this->k;
    xmin /= this->k;
    xmax /= this->k;

    return findValueBisectively(x, xmin, xmax);
}

double VirtualSineSurfaceFunction::findValueBisectively(double x, double xmin, double xmax) const {
    Expects(xmax > xmin);

    double currX{};
    double xmid{};
    do {
        xmid = (xmax + xmin) / 2;
        currX = calculateDiskCenterX(xmid);
        if (currX > x)
            xmax = xmid;
        else
            xmin = xmid;
    } while (xmax - xmin > PRECISION);

    return this->calculateDiskCenterY(xmid);
}

double VirtualSineSurfaceFunction::calculateDiskCenterX(double tangentX) const {
    return tangentX - A * std::cos(k*tangentX) * r * k / std::sqrt(1 + A*A*k*k*std::pow(std::cos(k*tangentX), 2));
}

double VirtualSineSurfaceFunction::calculateDiskCenterY(double tangentX) const {
    return A * std::sin(k * tangentX) + r / std::sqrt(1 + A * A * k * k * std::pow(std::cos(k*tangentX), 2));
}

VirtualSineSurfaceFunction::VirtualSineSurfaceFunction(double A, double k, double r) : A{A}, k{k}, r{r} {
    Expects(r > 0);
}

double VirtualSineSurfaceFunction::calculateMinValue() const {
    if (this->r < 1/this->A/this->k/this->k)
        return -this->A + this->r;

    double xmin = M_PI / this->k;
    double xmid = xmin;
    double xmax = 1.5 * M_PI / this->k;
    do {
        xmid = (xmid + xmax) / 2;
    } while (this->calculateDiskCenterX(xmid) < xmax);

    return this->findValueBisectively(xmax, xmin, xmid);
}
