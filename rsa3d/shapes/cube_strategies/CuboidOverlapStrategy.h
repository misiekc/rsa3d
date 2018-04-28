//
// Created by PKua on 03.12.17.
//

#ifndef RSA3D_CUBOIDOVERLAPSTRATEGY_H
#define RSA3D_CUBOIDOVERLAPSTRATEGY_H

#include "../Cuboid.h"

// Overlap stategies
class CuboidOverlapStrategy : public OverlapStrategy<3, 0> {
public:
    using VERTEX = Cuboid::VERTEX;
    using COORD = Cuboid::COORD;

    virtual std::string getName() const = 0;
    virtual void runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) const = 0;
};


#endif //RSA3D_OVERLAPSTRATEGY_H
