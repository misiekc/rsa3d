//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_OCTAHEDRON_H
#define RSA3D_OCTAHEDRON_H


#include "RegularSolid.h"
#include "../../OrderCalculable.h"

class Octahedron : public RegularSolid<Octahedron>, public OrderCalculable {
private:
    friend RegularSolid<Octahedron>;

    static void calculateStatic(const std::string &attr);

public:
    explicit Octahedron(const Matrix<3, 3> &orientation) : RegularSolid<Octahedron>{orientation} {};

    double projectionHalfsize(const Vector<3> &axis) const override;

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;
};


#endif //RSA3D_OCTAHEDRON_H
