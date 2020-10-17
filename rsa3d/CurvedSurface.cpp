//
// Created by Piotr Kubala on 17/10/2020.
//

#include "CurvedSurface.h"

void CurvedSurface::fillInLastCoordinate(RSAVector &position) const {
    this->surfaceFunction->fillInLastCoordinate(position);
}

SurfaceFunction::MinMax CurvedSurface::calculateValueRange(const RSAVector &voxelPosition,
                                                           double voxelSpatialSize) const
{
    return this->surfaceFunction->calculateValueRange(voxelPosition, voxelSpatialSize);
}
