//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARICOSAHEDRON_H
#define RSA3D_REGULARICOSAHEDRON_H


#include "PlatonicSolid.h"

class RegularIcosahedron : public PlatonicSolid<RegularIcosahedron> {
private:
    friend PlatonicSolid<RegularIcosahedron>;

    constexpr static double gold = (1 + std::sqrt(5.)) / 2;
    constexpr static double edge = std::pow(1.25 + 5. / 12 * std::sqrt(5.), -1./3);
    constexpr static double edgeFactor = edge / 2;

    constexpr static double circumsphereRadius = edge * std::sqrt(gold * std::sqrt(5)) / 2;
    constexpr static double insphereRadius = edge * gold * gold / 2 / std::sqrt(3);
    static void calculateStatic(const std::string &attr);

public:
    explicit RegularIcosahedron(const Matrix<3, 3> &orientation);

    double projectionHalfsize(const Vector<3> &axis) const;     /* CRTP implement */
    intersection::tri_polyh getTriangles() const;              /* CRTP implement */
};


#endif //RSA3D_REGULARICOSAHEDRON_H
