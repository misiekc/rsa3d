//
// Created by PKua on 07.02.18.
//

#ifndef RSA3D_SPHEROCYLINDER2D_H
#define RSA3D_SPHEROCYLINDER2D_H


#include "../AnisotropicShape2D.h"

class SpheroCylinder2D : public AnisotropicShape2D {
public:
    int pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) override;

    int pointInside(BoundaryConditions *bc, double *da) override;

    double getNeighbourListCellSize() override;

    double getVoxelSize() override;

    int overlap(BoundaryConditions *bc, Shape *s) override;

    double getVolume() override;

    std::string toString() override;

    std::string toPovray() const override;

    void store(std::ostream &f) const override;

    void restore(std::istream &f) override;

};


#endif //RSA3D_SPHEROCYLINDER2D_H
