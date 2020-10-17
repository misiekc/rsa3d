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
    std::unordered_set<int> activeGridCells;
    std::vector<int> randomAccessActiveGridCells;
    CurvedSurface *surface;

protected:
    bool isVoxelInsidePacking(Voxel *v) override;

public:
    CurvedSurfaceVoxelList(double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange,
                           double requestedAngularVoxelSize, CurvedSurface *surface);

    unsigned short splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl,
                               RSABoundaryConditions *bc) override;
    void getRandomEntry(RSAVector *position, RSAOrientation *orientation, Voxel **v, RND *rnd) override;
};


#endif //RSA3D_CURVEDSURFACEVOXELLIST_H
