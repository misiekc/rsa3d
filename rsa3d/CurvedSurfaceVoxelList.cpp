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

    this->activeGridCells.clear();
    this->activeGridCells.reserve(this->length);
    this->randomAccessActiveGridCells.clear();
    this->randomAccessActiveGridCells.reserve(this->length);

    for (std::size_t i{}; i < this->length; i++) {
        const auto &voxel = *(this->voxels[i]);
        int idx = this->getVoxelIndex(voxel);

        if (this->activeGridCells.find(idx) == this->activeGridCells.end()) {
            this->activeGridCells.insert(idx);
            this->randomAccessActiveGridCells.push_back(idx);
        }
    }

    return returnCode;
}

int CurvedSurfaceVoxelList::getVoxelIndex(const Voxel &voxel) const {
    const auto &pos = voxel.getPosition();
    double daArray[RSA_SPATIAL_DIMENSION];
    pos.copyToArray(daArray);
    double spatialSize = this->getSpatialVoxelSize();
    int n = (int)(spatialRange / spatialSize) + 1;
    int idx = position2i(daArray, RSA_SPATIAL_DIMENSION - 1, n*spatialSize, spatialSize, n);
    Assert(idx >= 0);
    return idx;
}

void CurvedSurfaceVoxelList::getRandomEntry(RSAVector *position, RSAOrientation *orientation, Voxel **v, RND *rnd) {
    if (!this->voxelsInitialized) {
        VoxelList::getRandomEntry(position, orientation, v, rnd);
        this->surface->fillInLastCoordinate(*position);
        Assert(std::abs((*position)[RSA_SPATIAL_DIMENSION - 1]) < this->spatialRange/2);
        (*position)[RSA_SPATIAL_DIMENSION - 1] += (this->spatialRange/2);
        return;
    }

    do {
        auto randomSetIdx = static_cast<std::size_t>(this->randomAccessActiveGridCells.size() * rnd->nextValue());
        Assert(randomSetIdx < this->randomAccessActiveGridCells.size());
        int randomGridCellIdx = this->randomAccessActiveGridCells[randomSetIdx];

        std::array<double, RSA_SPATIAL_DIMENSION> da{};
        da.fill(0);
        double spatialSize = this->getSpatialVoxelSize();
        int n = (int) (this->spatialRange / spatialSize) + 1;
        i2position(da.data(), RSA_SPATIAL_DIMENSION - 1, randomGridCellIdx, spatialSize, n);
        *position = da;

        bool insidePacking = true;
        for (std::size_t i{}; i < RSA_SPATIAL_DIMENSION - 1; i++) {
            (*position)[i] += ((rnd->nextValue() - 0.5) * spatialSize);
            if ((*position)[i] >= this->spatialRange) {
                insidePacking = false;
                break;
            }
        }
        if (!insidePacking)
            continue;

        this->surface->fillInLastCoordinate(*position);
        Assert(std::abs((*position)[RSA_SPATIAL_DIMENSION - 1]) < this->spatialRange / 2);
        (*position)[RSA_SPATIAL_DIMENSION - 1] += (this->spatialRange / 2);

        RSAOrientation dummyOrientation{};
        (*v) = this->getVoxel(*position, dummyOrientation);
        if (*v == nullptr) {
            std::cout << "no voxel - setIdx: " << randomSetIdx << ", gridIdx: " << randomGridCellIdx;
            std::cout << ", position: " << *position << ", spatial size: " << spatialSize << ", candidates: " << std::endl;
            bool rangeShown = false;
            for (std::size_t i{}; i < this->length; i++) {
                const auto &voxel = *(this->voxels[i]);
                if (this->getVoxelIndex(voxel) == randomGridCellIdx) {
                    if (!rangeShown) {
                        auto [min, max] = this->surface->calculateValueRange(voxel.getPosition(), spatialSize);
                        min += this->spatialRange/2;
                        max += this->spatialRange/2;
                        std::cout << "  function range: {" << min << ", " << max << "}" << std::endl;
                        rangeShown = true;
                    };
                    std::cout << "  {" << voxel.getPosition() << ", " << (voxel.getPosition() + RSAVector(spatialSize));
                    std::cout << "}";

                    if (this->isVoxelInsidePacking(&voxel)) {
                        std::cout << std::endl;
                    } else {
                        std::cout << " - not inside packing";
                        auto it = this->registeredVoxels.find(&voxel);
                        if (it == this->registeredVoxels.end()) {
                            std::cout << " - not registered" << std::endl;
                        } else {
                            std::cout << " - registered as: " << std::endl << "    ";
                            std::copy(it->second.begin(), it->second.end(), std::ostream_iterator<VoxelEntry>(std::cout, "\n    "));
                            std::cout << std::endl;
                        }
                    }
                }
            }
        }
    } while (*v == nullptr);
}

bool CurvedSurfaceVoxelList::isVoxelInsidePacking(const Voxel *v) const {
    this->registeredVoxels[v].push_back({v->getPosition(), this->getSpatialVoxelSize()});

    if (!VoxelList::isVoxelInsidePacking(v))
        return false;

    RSAVector voxelPos = v->getPosition();
    double spatialSize = this->getSpatialVoxelSize();
    auto [min, max] = this->surface->calculateValueRange(voxelPos, spatialSize);
    min += this->spatialRange/2;
    max += this->spatialRange/2;
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
