//
// Created by Piotr Kubala on 17/10/2020.
//

#ifndef RSA3D_CURVEDSURFACEVOXELLIST_H
#define RSA3D_CURVEDSURFACEVOXELLIST_H

#include <unordered_set>
#include <memory>

#include "VoxelList.h"
#include "CurvedSurface.h"


class CurvedSurfaceVoxelList : public VoxelList {
private:
    struct VoxelEntry {
        RSAVector position;
        double size{};

        friend std::ostream &operator<<(std::ostream &out, const VoxelEntry &entry) {
            return out << "{" << entry.position << ", " << entry.size << "}";
        }
    };

    std::unordered_set<int> activeGridCells;
    std::vector<int> randomAccessActiveGridCells;
    CurvedSurface *surface;
    mutable std::unordered_map<const Voxel *, std::vector<VoxelEntry>> registeredVoxels;

    int getVoxelIndex(const Voxel &voxel) const;

protected:
    bool isVoxelInsidePacking(const Voxel *v) const override;

public:
    CurvedSurfaceVoxelList(double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange,
                           double requestedAngularVoxelSize, CurvedSurface *surface);

    unsigned short splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl,
                               RSABoundaryConditions *bc) override;
    void getRandomEntry(RSAVector *position, RSAOrientation *orientation, Voxel **v, RND *rnd) override;
};


#endif //RSA3D_CURVEDSURFACEVOXELLIST_H
