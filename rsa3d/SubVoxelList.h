//
// Created by Michal Ciesla on 05/03/2026.
//

#ifndef RSA3D_TOPLEVELVOXELVOXELLIST_H
#define RSA3D_TOPLEVELVOXELVOXELLIST_H

#include "VoxelList.h"

class SubVoxelList: public VoxelList{

private:
    size_t initialLength;
    size_t *parents;
    bool *nulls;
    void testList();

protected:
    void moveVoxelInList(size_t from, size_t to) override;

public:
    ~SubVoxelList() override;
    SubVoxelList(const VoxelList &vl, size_t beginIndex, size_t endIndex);
    SubVoxelList(const VoxelList &vl, size_t topIndex);

    unsigned short splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc, bool verbose=true) override;

    [[nodiscard]] std::vector<size_t> getIndicesOfRemovedVoxels() const;



};


#endif //RSA3D_TOPLEVELVOXELVOXELLIST_H