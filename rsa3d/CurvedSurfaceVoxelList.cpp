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
        const auto &voxel = *(this->voxels[i]);
        int idx = getSurfaceCellIndexForVoxel(voxel);

        if (this->activeSurfaceCells.find(idx) == this->activeSurfaceCells.end()) {
            this->activeSurfaceCells.insert(idx);
            this->randomAccessActiveSurfaceCells.push_back(idx);
        }
    }
}

int CurvedSurfaceVoxelList::getSurfaceCellIndexForVoxel(const Voxel &voxel) const {
    const auto &pos = voxel.getPosition();
    double posArray[RSA_SPATIAL_DIMENSION];
    pos.copyToArray(posArray);
    double spatialSize = this->getSpatialVoxelSize();
    int n = (int)(spatialRange / spatialSize) + 1;
    int idx = position2i(posArray, RSA_SPATIAL_DIMENSION - 1, n * spatialSize, spatialSize, n);
    Assert(idx >= 0);
    return idx;
}

void CurvedSurfaceVoxelList::getRandomEntry(RSAVector *position, RSAOrientation *orientation, Voxel **v, RND *rnd) {
    if (!this->voxelsInitialized) {
        VoxelList::getRandomEntry(position, orientation, v, rnd);
        this->fillInLastCoordinate(*position);
        return;
    }

    do {
        // Sample random active grid cell
        auto randomSetIdx = static_cast<std::size_t>(this->randomAccessActiveSurfaceCells.size() * rnd->nextValue());
        Assert(randomSetIdx < this->randomAccessActiveSurfaceCells.size());
        int randomGridCellIdx = this->randomAccessActiveSurfaceCells[randomSetIdx];

        // Calculate position of its CENTER (as i2position does) and store in *position
        *position = this->calculateSurfaceCellBottomLeftPosition(randomGridCellIdx);

        // Sample random position inside this cell and store in *position
        for (std::size_t i{}; i < RSA_SPATIAL_DIMENSION - 1; i++)
            (*position)[i] += (rnd->nextValue() * this->getSpatialVoxelSize());

        // Fill in last coordinate according to surface function
        this->fillInLastCoordinate(*position);

        // Find voxel contatining this point
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

    // Out surface is in the middle
    min += this->spatialRange/2;
    max += this->spatialRange/2;

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

    // Surface places a function along last coordinate = 0; We want it to be in the middle of the simulation box
    Assert(std::abs(position[RSA_SPATIAL_DIMENSION - 1]) < this->spatialRange/2);
    position[RSA_SPATIAL_DIMENSION - 1] += (this->spatialRange/2);
}

double CurvedSurfaceVoxelList::getVoxelsVolume() {
    if (!this->voxelsInitialized)
        return std::pow(this->spatialRange, RSA_SPATIAL_DIMENSION - 1);

    // The volume is only calculated among all but the last one axis, being the function value axis, so we are using
    // the surface grid, not the voxels.
    // This is because points are sampled in all directions apart from that last one.
    double result = 0;
    double spatialSize = this->getSpatialVoxelSize();
    for (int cellIdx : randomAccessActiveSurfaceCells) {
        double s = 1.0;
        RSAVector position = this->calculateSurfaceCellBottomLeftPosition(cellIdx);
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

RSAVector CurvedSurfaceVoxelList::calculateSurfaceCellBottomLeftPosition(int surfaceCellIdx) {
    RSAVector result;
    std::array<double, RSA_SPATIAL_DIMENSION> positionArray{};
    positionArray.fill(0);
    double spatialSize = this->getSpatialVoxelSize();
    int n = (int) (this->spatialRange / spatialSize) + 1;
    i2position(positionArray.data(), RSA_SPATIAL_DIMENSION - 1, surfaceCellIdx, spatialSize, n);
    result = positionArray;
    result -= RSAVector(spatialSize / 2);
    return result;
}
