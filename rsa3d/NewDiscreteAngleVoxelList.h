//
// Created by Michal Ciesla on 12/23/22.
//

#ifndef RSA3D_NEWDISCRETEANGLEVOXELLIST_H
#define RSA3D_NEWDISCRETEANGLEVOXELLIST_H


#include "VoxelList.h"

#include <unordered_map>

class NewDiscreteAngleVoxelList : public VoxelList{
private:
    std::uniform_real_distribution<double> *u01Distribution;

    double *voxelMap;

    std::vector<Orientation<1>> allowedOrientations;
    std::vector<Orientation<1>> *orientationsVector;
    Orientation<1> orientationsVectorRange;
    std::vector<Orientation<1>> *orientationsHalfVector;
    double orientationsHalfVectorRange;

    Orientation<1> getAngularVoxelSize(size_t i) const;
    double getVoxelVolume(size_t i) const;
    std::vector<Orientation<1>> getPossibleOrientations(double start, double range) const;
    void createVoxelMap();
    void createOrientationsMap(double range);

protected:
    virtual void allocateVoxels(size_t size) override;

public:
    NewDiscreteAngleVoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, const Orientation<1> &shapeAngularRange, const Orientation<1> &requestedAngularVoxelSize, const std::vector<Orientation<1>> &orientations);
    virtual ~NewDiscreteAngleVoxelList();

    void getRandomEntry(RSAVector *position, Orientation<1> *orientation, Voxel **v, RND *rnd);
    void getRandomPositionAndOrientation(RSAVector *position, Orientation<1> *orientation, Voxel *v, RND *rnd);

    unsigned short splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc) override;
    double getVoxelsVolume() const override;

    bool analyzeVoxel(Voxel *v, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc, double spatialSize, const Orientation<1> &angularSize, unsigned short depth=0) const override;

    };
#endif //RSA3D_NEWDISCRETEANGLEVOXELLIST_H
