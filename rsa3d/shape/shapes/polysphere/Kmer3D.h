//
// Created by ciesla on 19.11.2025.
//

#ifndef RSA3D_KMER3D_H
#define RSA3D_KMER3D_H


#include "Polysphere.h"

class Kmer3D : public Polysphere {
protected:
    std::array<Matrix<3, 3>, 2> getMinMaxMatrices(const Orientation<3> &voxelOrientation, const Orientation<3> &angularSize) const override;

public:
    static void initClass(const std::string &attr);

    static double calculateVolume(int numberOfDisks, double distance);
    static std::string preparePolysphereAttr(int numberOfDisks, double length);
    bool voxelInside(BoundaryConditions<3> *bc, const Vector<3> &voxelPosition, const Orientation<3> &voxelOrientation, double spatialSize, const Orientation<3> &angularSize) const override;
};

#endif //RSA3D_KMER_H
