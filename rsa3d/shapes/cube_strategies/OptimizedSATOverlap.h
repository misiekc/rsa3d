//
// Created by PKua on 05.12.17.
//

#ifndef RSA3D_OPTIMIZEDSATOVERLAP_H
#define RSA3D_OPTIMIZEDSATOVERLAP_H


#include "OverlapStrategy.h"

class OptimizedSATOverlap : public OverlapStrategy {
public:
    bool overlap(const Cuboid *cube1, const Cuboid *cube2, BoundaryConditions *bc) override;

    std::string getName() override;

    void runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) override;
};


#endif //RSA3D_OPTIMIZEDSATOVERLAP_H
