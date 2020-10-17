//
// Created by Piotr Kubala on 17/10/2020.
//

#include "CurvedSurfaceVoxelList.h"
#include "utils/Utils.h"

unsigned short CurvedSurfaceVoxelList::splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl,
                                                   RSABoundaryConditions *bc)
{
    auto returnCode = VoxelList::splitVoxels(minDx, maxVoxels, nl, bc);

    this->activeGridCells.clear();
    this->activeGridCells.reserve(this->length);
    this->randomAccessActiveGridCells.clear();
    this->randomAccessActiveGridCells.reserve(this->length);

    for (std::size_t i{}; i < this->length; i++) {
        const auto &voxel = *(this->voxels[i]);
        const auto &pos = voxel.getPosition();
        double daArray[RSA_SPATIAL_DIMENSION];
        pos.copyToArray(daArray);
        double spatialSize = this->getSpatialVoxelSize();
        int n = (int)(this->spatialRange/spatialSize) + 1;
        int idx = position2i(daArray, RSA_SPATIAL_DIMENSION, n*spatialSize, spatialSize, n);

        if (this->activeGridCells.find(idx) != this->activeGridCells.end()) {
            this->activeGridCells.insert(idx);
            this->randomAccessActiveGridCells.push_back(idx);
        }
    }

    return returnCode;
}

void CurvedSurfaceVoxelList::getRandomEntry(RSAVector *position, RSAOrientation *orientation, Voxel **v, RND *rnd) {
    if (!this->voxelsInitialized) {
        VoxelList::getRandomEntry(position, orientation, v, rnd);
        return;
    }

    auto randomSetIdx = static_cast<std::size_t>(this->randomAccessActiveGridCells.size() * rnd->nextValue());
    Assert(randomSetIdx < this->randomAccessActiveGridCells.size());
    int randomGridCellIdx = this->randomAccessActiveGridCells[randomSetIdx];

    std::array<double, RSA_SPATIAL_DIMENSION> da{};
    da.fill(0);
    double spatialSize = this->getSpatialVoxelSize();
    int n = (int)(this->spatialRange/spatialSize) + 1;
    i2position(da.data(), RSA_SPATIAL_DIMENSION, randomGridCellIdx, spatialSize, n);
    *position = da;
    for (std::size_t i{}; i < RSA_SPATIAL_DIMENSION - 1; i++)
        (*position)[i] += (rnd->nextValue() * spatialSize);
    this->surface->fillInLastCoordinate(*position);
    Assert(std::abs((*position)[RSA_SPATIAL_DIMENSION - 1]) < this->spatialRange/2);
    (*position)[RSA_SPATIAL_DIMENSION - 1] += (this->spatialRange/2);

    RSAOrientation dummyOrientation{};
    (*v) = this->getVoxel(*position, dummyOrientation);
}

bool CurvedSurfaceVoxelList::isVoxelInsidePacking(Voxel *v) {
    if (!VoxelList::isVoxelInsidePacking(v))
        return false;

    RSAVector voxelPos = v->getPosition();
    double spatialSize = this->getSpatialVoxelSize();
    auto [min, max] = this->surface->calculateValueRange(voxelPos, spatialSize);
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
    Expects(RSA_ANGULAR_DIMENSION == 0);
    Expects(RSA_SPATIAL_DIMENSION > 1);
}
