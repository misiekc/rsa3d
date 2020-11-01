//
// Created by Piotr Kubala on 17/10/2020.
//

#include <iterator>

#include "CurvedSurfaceVoxelList.h"
#include "utils/Utils.h"

unsigned short CurvedSurfaceVoxelList::splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl,
                                                   RSABoundaryConditions *bc)
{
    auto returnCode = VoxelList::splitVoxels(minDx, maxVoxels, nl, bc);
    this->rebuildActiveSurfaceCells();
    return returnCode;
}

void CurvedSurfaceVoxelList::rebuildActiveSurfaceCells() {
    this->activeSurfaceCells.clear();
    this->activeSurfaceCells.reserve(this->length);
    this->randomAccessActiveSurfaceCells.clear();
    this->randomAccessActiveSurfaceCells.reserve(this->length);

    for (std::size_t i{}; i < this->length; i++) {
        RSAVector pos = this->voxels[i]->getPosition();
        pos[RSA_SPATIAL_DIMENSION - 1] = 0; // 2 voxels differing only on last coordinate should be considered the same

        if (this->activeSurfaceCells.find(pos) == this->activeSurfaceCells.end()) {
            this->activeSurfaceCells.insert(pos);
            this->randomAccessActiveSurfaceCells.push_back(pos);
        }
    }
}

void CurvedSurfaceVoxelList::getRandomEntry(RSAVector *position, RSAOrientation *orientation, Voxel **v, RND *rnd) {
    if (!this->voxelsInitialized) {
        VoxelList::getRandomEntry(position, orientation, v, rnd);
        this->fillInLastCoordinate(*position);
        return;
    }

    do {
        auto randomSetIdx = static_cast<std::size_t>(this->randomAccessActiveSurfaceCells.size() * rnd->nextValue());
        Assert(randomSetIdx < this->randomAccessActiveSurfaceCells.size());
        *position = this->randomAccessActiveSurfaceCells[randomSetIdx];

        // Sample random position inside this cell and store in *position
        for (std::size_t i{}; i < RSA_SPATIAL_DIMENSION - 1; i++)
            (*position)[i] += (rnd->nextValue() * this->getSpatialVoxelSize());

        this->fillInLastCoordinate(*position);
        RSAOrientation dummyOrientation{};
        (*v) = this->getVoxel(*position, dummyOrientation);

        // Repeat until a voxel is found - if a grid cell is active that means there are still voxels in it, so sooner
        // or later the voxel will be found (unless is is outside of the surface and not deleted...)
    } while (*v == nullptr);
}

bool CurvedSurfaceVoxelList::isVoxelInsidePacking(const Voxel *v, double spatialSize) const {
    Expects(spatialSize > 0);

    if (!VoxelList::isVoxelInsidePacking(v, spatialSize))
        return false;

    RSAVector voxelPos = v->getPosition();
    auto [min, max] = this->surface->calculateValueRange(voxelPos, spatialSize);

    // We place the surface in last coordinate = 0
    min -= this->surface->getValueSpan().min;
    max -= this->surface->getValueSpan().min;

    // If voxel is below or above the surface, it should be deleted
    if (max < voxelPos[RSA_SPATIAL_DIMENSION - 1] || min > voxelPos[RSA_SPATIAL_DIMENSION - 1] + spatialSize)
        return false;
    else
        return true;
}

CurvedSurfaceVoxelList::CurvedSurfaceVoxelList(double packingSpatialSize, double requestedSpatialVoxelSize,
                                               double shapeAngularRange, double requestedAngularVoxelSize,
                                               CurvedSurface *surface)
        : VoxelList(RSA_SPATIAL_DIMENSION, packingSpatialSize, requestedSpatialVoxelSize, shapeAngularRange,
                    requestedAngularVoxelSize),
          surface(surface)
{
    // Currently we are supporting only isotropic shapes
    Expects(RSA_ANGULAR_DIMENSION == 0);
    // Function plot is at least 2-dim
    Expects(RSA_SPATIAL_DIMENSION > 1);
}

void CurvedSurfaceVoxelList::fillInLastCoordinate(RSAVector &position) {
    this->surface->fillInLastCoordinate(position);

    // Minimal value of function should land on 0 on the last coordinate
    position[RSA_SPATIAL_DIMENSION - 1] -= this->surface->getValueSpan().min;
    Ensures(position[RSA_SPATIAL_DIMENSION - 1] >= 0);
}

double CurvedSurfaceVoxelList::getVoxelsVolume() {
    if (!this->voxelsInitialized)
        return std::pow(this->spatialRange, RSA_SPATIAL_DIMENSION - 1);

    // The volume is only calculated among all but the last one axis, being the function value axis, so we are using
    // the surface grid, not the voxels.
    // This is because points are sampled in all directions apart from that last one.
    double result = 0;
    double spatialSize = this->getSpatialVoxelSize();
    for (const auto &position : randomAccessActiveSurfaceCells) {
        double s = 1.0;
        for (std::size_t i{}; i < RSA_SPATIAL_DIMENSION - 1; i++) {
            if (position[i] + spatialSize > this->spatialRange)
                s *= this->spatialRange - position[i];
            else
                s *= spatialSize;
        }
        result += s;
    }
    return result;
}

std::array<std::size_t, RSA_SPATIAL_DIMENSION> CurvedSurfaceVoxelList::calculateVoxelsGridLinearSize() const {
    auto result = VoxelList::calculateVoxelsGridLinearSize();
    double span = this->surface->getValueSpan().getSpan();
    result[RSA_SPATIAL_DIMENSION - 1] = this->findArraySize(span, this->initialVoxelSize);
    return result;
}
