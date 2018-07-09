//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARDODECAHEDRON_H
#define RSA3D_REGULARDODECAHEDRON_H


#include "PlatonicSolid.h"
#include "../../OrderCalculable.h"

class RegularDodecahedron : public PlatonicSolid<RegularDodecahedron>, public OrderCalculable {
private:
    friend PlatonicSolid<RegularDodecahedron>;

    constexpr static double gold = (1 + std::sqrt(5.)) / 2;
    constexpr static double edge = std::pow(3.75 + 1.75 * std::sqrt(5.), -1./3);
    constexpr static double edgeFactor = edge * gold / 2;

    constexpr static double circumsphereRadius = edge * std::sqrt(3.) / 2 * gold;
    constexpr static double insphereRadius = edge * gold * gold / 2 / std::sqrt(3 - gold);

    static void calculateStatic(const std::string &attr);

public:
    explicit RegularDodecahedron(const Matrix<3, 3> &orientation);

    double projectionHalfsize(const Vector<3> &axis) const;     /* CRTP implement */
    intersection::tri_polyh getTriangles() const;              /* CRTP implement */

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;
};


#endif //RSA3D_REGULARDODECAHEDRON_H
