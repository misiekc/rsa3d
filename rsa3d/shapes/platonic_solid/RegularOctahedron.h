//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULAROCTAHEDRON_H
#define RSA3D_REGULAROCTAHEDRON_H


#include "PlatonicSolid.h"

class RegularOctahedron : public PlatonicSolid<RegularOctahedron> {
private:
    friend PlatonicSolid<RegularOctahedron>;

    constexpr static double exsphereRadius = 0;
    constexpr static double insphereRadius = 0;
    static std::array<Vector<3>, 4> faceAxes;
    static std::array<Vector<3>, 6> edgeAxes;

    static void calculateStatic(const std::string &attr) {

    }

public:
    explicit RegularOctahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularOctahedron>{orientation} {}

    double projectionHalfsize(const Vector<3> &axis) const {
        return 0;
    }
};


#endif //RSA3D_REGULAROCTAHEDRON_H
