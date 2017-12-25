//
// Created by PKua on 03.12.17.
//

#ifndef RSA3D_MINEOVERLAP_H
#define RSA3D_MINEOVERLAP_H

#include "OverlapStrategy.h"

class MineOverlap : public OverlapStrategy {
private:
    bool checkSegment(Cuboid *cube, const Vector<3> & point1, const Vector<3> & point2);

public:
    bool overlap(Cuboid *cube1, Cuboid *cube2, BoundaryConditions *bc) override;
    std::string getName() override;

    void runOverheadOperations(Cuboid *cube1, Cuboid *cube2) override;
};

#endif //RSA3D_MINEOVERLAP_H
