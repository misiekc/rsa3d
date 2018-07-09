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

    static void calculateStatic(const std::string &attr);

public:
    explicit RegularDodecahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularDodecahedron>{orientation} {};

    double projectionHalfsize(const Vector<3> &axis) const;     /* CRTP implement */

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;
};


#endif //RSA3D_REGULARDODECAHEDRON_H
