//
// Created by PKua on 03.12.17.
//

#ifndef RSA3D_TRITRIOVERLAP_H
#define RSA3D_TRITRIOVERLAP_H

#include "OverlapStrategy.h"

class TriTriOverlap : public OverlapStrategy {
private:
    void obtainTris(Cuboid * cube, Vector<3> (&arr)[12][3], const Vector<3> & translation);

public:
    bool overlap(Cuboid *cube1, Cuboid *cube2, BoundaryConditions *bc) override;
    std::string getName() override;

    void runOverheadOperations(Cuboid *cube1, Cuboid *cube2) override;
};


#endif //RSA3D_TRITRIOVERLAP_H
