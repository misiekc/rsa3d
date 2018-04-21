//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARICOSAHEDRON_H
#define RSA3D_REGULARICOSAHEDRON_H


#include "PlatonicSolid.h"

class RegularIcosahedron : public PlatonicSolid<RegularIcosahedron> {
private:
    friend PlatonicSolid<RegularIcosahedron>;

    constexpr static double exsphereRadius = 0;
    constexpr static double insphereRadius = 0;

    static void calculateStatic(const std::string &attr) {

    }

public:
    explicit RegularIcosahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularIcosahedron>{orientation} {}
};


#endif //RSA3D_REGULARICOSAHEDRON_H
