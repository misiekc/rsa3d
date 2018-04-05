//
// Created by PKua on 03.12.17.
//

#ifndef RSA3D_TRITRIOVERLAP_H
#define RSA3D_TRITRIOVERLAP_H

#include "OverlapStrategy.h"

class TriTriOverlap : public OverlapStrategy {
private:
    void obtainTris(const Cuboid *cube, Vector<3> (&arr)[12][3], const Vector<3> &translation);

public:
    bool overlap(const Cuboid *cube1, const Cuboid *cube2, BoundaryConditions *bc) override;
    std::string getName() override;

    void runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) override;
};


#endif //RSA3D_TRITRIOVERLAP_H
