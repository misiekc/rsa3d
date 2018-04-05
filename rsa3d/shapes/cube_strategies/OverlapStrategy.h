//
// Created by PKua on 03.12.17.
//

#ifndef RSA3D_OVERLAPSTRATEGY_H
#define RSA3D_OVERLAPSTRATEGY_H

#include "../Cuboid.h"

// Overlap stategies
class OverlapStrategy {
public:
    using VERTEX = Cuboid::VERTEX;
    using COORD = Cuboid::COORD;

    virtual ~OverlapStrategy() = default;
    virtual bool overlap(const Cuboid * cube1, const Cuboid * cube2, BoundaryConditions *bc) = 0;
    virtual std::string getName() = 0;
    virtual void runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) = 0;
};


#endif //RSA3D_OVERLAPSTRATEGY_H
