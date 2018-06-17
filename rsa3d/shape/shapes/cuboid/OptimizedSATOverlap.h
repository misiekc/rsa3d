//
// Created by PKua on 05.12.17.
//

#ifndef RSA3D_OPTIMIZEDSATOVERLAP_H
#define RSA3D_OPTIMIZEDSATOVERLAP_H


#include "CuboidOverlapStrategy.h"

class OptimizedSATOverlap : public CuboidOverlapStrategy {
public:
    bool overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const override;
    std::string getName() const override;
    void runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) const override;
};


#endif //RSA3D_OPTIMIZEDSATOVERLAP_H
