//
// Created by Piotr Kubala on 17/10/2020.
//

#ifndef RSA3D_CURVEDSURFACEVOXELLIST_H
#define RSA3D_CURVEDSURFACEVOXELLIST_H

#include <unordered_set>
#include <memory>

#include "VoxelList.h"
#include "CurvedSurface.h"

/**
 * @brief A VoxelList for CurvedSurface.
 * @details It samples coordinates randomly in all but last coordinates and the last coordinate is given by a
 * SurfaceFunction value in this point. Also, it removes the voxels which do not overlap the surface.
 */
class CurvedSurfaceVoxelList : public VoxelList {
private:
    // Debug helper struct used by registeredVoxels
    struct VoxelEntry {
        RSAVector position;
        double size{};

        friend std::ostream &operator<<(std::ostream &out, const VoxelEntry &entry) {
            return out << "{" << entry.position << ", " << entry.size << "}";
        }
    };

    // Those 2 structures are used to track which parts of surface have active voxels. They contain indices computed
    // by position2i for a grid like voxels with current voxel spatial size, but without last coordinate.
    // So if a cell representing 2d position {1.5, 3} is active (for spatial size 0.5), it means that for example there
    // are 2 active voxels with coordinates {1.5, 3, 4} and {1.4, 3, 4.5}
    std::unordered_set<int> activeSurfaceCells;
    std::vector<int> randomAccessActiveSurfaceCells;

    CurvedSurface *surface;

    // A debug map registering voxels comming to isVoxelInsidePacking
    mutable std::unordered_map<const Voxel *, std::vector<VoxelEntry>> registeredVoxels;

    int getSurfaceCellIndexForVoxel(const Voxel &voxel) const;
    void rebuildActiveSurfaceCells();
    void fillInLastCoordinate(RSAVector &position);

protected:
    bool isVoxelInsidePacking(const Voxel *v, double spatialSize) const override;

public:
    CurvedSurfaceVoxelList(double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange,
                           double requestedAngularVoxelSize, CurvedSurface *surface);

    unsigned short splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl,
                               RSABoundaryConditions *bc) override;
    void getRandomEntry(RSAVector *position, RSAOrientation *orientation, Voxel **v, RND *rnd) override;
};


#endif //RSA3D_CURVEDSURFACEVOXELLIST_H
