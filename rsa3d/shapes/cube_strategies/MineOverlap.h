//
// Created by PKua on 03.12.17.
//

#ifndef RSA3D_MINEOVERLAP_H
#define RSA3D_MINEOVERLAP_H

#include "CuboidOverlapStrategy.h"

class MineOverlap : public CuboidOverlapStrategy {
private:
    bool checkSegment(const Cuboid *cube, const Vector<3> &point1, const Vector<3> &point2);

public:
    bool overlap(const Cuboid *cube1, const Cuboid *cube2, BoundaryConditions *bc) override;
    std::string getName() override;

    void runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) override;
};

#endif //RSA3D_MINEOVERLAP_H
