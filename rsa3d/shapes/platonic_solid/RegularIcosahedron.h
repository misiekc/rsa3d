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
    static std::array<Vector<3>, 10> faceAxes;
    static std::array<Vector<3>, 15> edgeAxes;

    static void calculateStatic(const std::string &attr) {

    }

public:
    explicit RegularIcosahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularIcosahedron>{orientation} {}

    double projectionHalfsize(const Vector<3> &axis) const {
        return 0;
    }
};


#endif //RSA3D_REGULARICOSAHEDRON_H
