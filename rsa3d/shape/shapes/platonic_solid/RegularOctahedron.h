//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULAROCTAHEDRON_H
#define RSA3D_REGULAROCTAHEDRON_H


#include "PlatonicSolid.h"

class RegularOctahedron : public PlatonicSolid<RegularOctahedron> {
private:
    friend PlatonicSolid<RegularOctahedron>;

    constexpr static double circumsphereRadius = std::pow(0.75, 1./3);
    constexpr static double insphereRadius = std::pow(48, -1./6);
    constexpr static double edgeFactor = std::pow(0.75, 1./3);

    static void calculateStatic(const std::string &attr);

public:
    explicit RegularOctahedron(const Matrix<3, 3> &orientation);

    double projectionHalfsize(const Vector<3> &axis) const;     /* CRTP implement */
    intersection::tri_polyh getTriangles() const;               /* CRTP implement */
};


#endif //RSA3D_REGULAROCTAHEDRON_H
