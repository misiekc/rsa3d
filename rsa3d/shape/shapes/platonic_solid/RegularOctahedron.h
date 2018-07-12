//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULAROCTAHEDRON_H
#define RSA3D_REGULAROCTAHEDRON_H


#include "PlatonicSolid.h"
#include "../../OrderCalculable.h"

class RegularOctahedron : public PlatonicSolid<RegularOctahedron>, public OrderCalculable {
private:
    friend PlatonicSolid<RegularOctahedron>;

    static void calculateStatic(const std::string &attr);

public:
    explicit RegularOctahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularOctahedron>{orientation} {};

    double projectionHalfsize(const Vector<3> &axis) const;         /* CRTP implement */

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;
};


#endif //RSA3D_REGULAROCTAHEDRON_H
