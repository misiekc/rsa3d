//
// Created by Piotr Kubala on 17/10/2020.
//

#include "CurvedSurface.h"

void CurvedSurface::fillInLastCoordinate(RSAVector &position) const {
    this->surfaceFunction->fillInLastCoordinate(position);
    Ensures(position[RSA_SPATIAL_DIMENSION - 1] >= this->valueSpan.min);
    Ensures(position[RSA_SPATIAL_DIMENSION - 1] <= this->valueSpan.max);
}

SurfaceFunction::MinMax CurvedSurface::calculateValueRange(const RSAVector &voxelPosition,
                                                           double voxelSpatialSize) const
{
    return this->surfaceFunction->calculateValueRange(voxelPosition, voxelSpatialSize);
}

double CurvedSurface::getArea() const {
    return std::pow(this->size, RSA_SPATIAL_DIMENSION - 1);
}

SurfaceFunction::MinMax CurvedSurface::getValueSpan() const {
    return this->valueSpan;
}

CurvedSurface::CurvedSurface(int dim, double s, double ndx, double vdx,
                             std::unique_ptr<SurfaceFunction> surfaceFunction,
                             std::unique_ptr<RSABoundaryConditions> bc)
        : Surface(dim, s, ndx, vdx, std::move(bc)), surfaceFunction(std::move(surfaceFunction))
{
    this->valueSpan = this->surfaceFunction->calculateValueRange(RSAVector{}, s);
    Expects(this->valueSpan.getSpan() < s);
}
