//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARDODECAHEDRON_H
#define RSA3D_REGULARDODECAHEDRON_H


#include "PlatonicSolid.h"

class RegularDodecahedron : public PlatonicSolid<RegularDodecahedron> {
private:
    friend PlatonicSolid<RegularDodecahedron>;

    constexpr static double exsphereRadius = 0;
    constexpr static double insphereRadius = 0;
    static std::array<Vector<3>, 6> faceAxes;
    static std::array<Vector<3>, 15> edgeAxes;

    static void calculateStatic(const std::string &attr) {

    }

public:
    explicit RegularDodecahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularDodecahedron>{orientation} {}

    double projectionHalfsize(const Vector<3> &axis) const {
        return 0;
    }
};


#endif //RSA3D_REGULARDODECAHEDRON_H
