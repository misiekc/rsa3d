//
// Created by PKua on 05.12.17.
//

#ifndef RSA3D_OPTIMIZEDSATOVERLAP_H
#define RSA3D_OPTIMIZEDSATOVERLAP_H


#include "OverlapStrategy.h"

class OptimizedSATOverlap : public OverlapStrategy {
public:
    bool overlap(Cuboid *cube1, Cuboid *cube2, BoundaryConditions *bc) override;

    std::string getName() override;
};


#endif //RSA3D_OPTIMIZEDSATOVERLAP_H
