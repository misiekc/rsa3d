//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARTETRAHEDRON_H
#define RSA3D_REGULARTETRAHEDRON_H


#include "PlatonicSolid.h"

class RegularTetrahedron : public PlatonicSolid<RegularTetrahedron> {
private:
    friend PlatonicSolid<RegularTetrahedron>;

    constexpr static double exsphereRadius = 0;
    constexpr static double insphereRadius = 0;

    static void calculateStatic(const std::string &attr) {

    }

public:
    explicit RegularTetrahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularTetrahedron>{orientation} {}
};


#endif //RSA3D_REGULARTETRAHEDRON_H
