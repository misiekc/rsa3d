//
// Created by PKua on 03.12.17.
//

#ifndef RSA3D_TRITRIOVERLAP_H
#define RSA3D_TRITRIOVERLAP_H

#include "CuboidOverlapStrategy.h"
#include "../../Intersection.h"

class TriTriOverlap : public CuboidOverlapStrategy {
public:
    int overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const override;
    std::string getName() const override;

    void runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) const override;
};


#endif //RSA3D_TRITRIOVERLAP_H
