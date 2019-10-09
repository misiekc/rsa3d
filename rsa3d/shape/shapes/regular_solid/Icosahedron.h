//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_ICOSAHEDRON_H
#define RSA3D_ICOSAHEDRON_H


#include "RegularSolid.h"

class Icosahedron : public RegularSolid<Icosahedron> {
private:
    friend RegularSolid<Icosahedron>;

    static ShapeData calculateStatic(const std::string &attr);

public:
    explicit Icosahedron(const Matrix<3, 3> &orientation) : RegularSolid<Icosahedron>{orientation} {};

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;

    double projectionHalfsize(const Vector<3> &axis) const override;
};


#endif //RSA3D_ICOSAHEDRON_H
