//
// Created by PKua on 12.07.18.
//

#ifndef RSA3D_SNUBCUBE_H
#define RSA3D_SNUBCUBE_H


#include "RegularSolid.h"
#include "UnoptimizedSATOverlapRS.h"
#include "Octahedron.h"

class SnubCube : public RegularSolid<SnubCube> {
private:
    friend RegularSolid<SnubCube>;

    using SymmetryPlatonicSolid = Octahedron;

    const static UnoptimizedSATOverlapRS overlapStrategy;

    static ShapeData calculateStatic(const std::string &attr);

public:
    explicit SnubCube(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;
};


#endif //RSA3D_SNUBCUBE_H
