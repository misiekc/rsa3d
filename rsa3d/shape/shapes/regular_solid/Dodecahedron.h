//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_DODECAHEDRON_H
#define RSA3D_DODECAHEDRON_H


#include "RegularSolid.h"
#include "../../OrderCalculable.h"

class Dodecahedron : public RegularSolid<Dodecahedron>, public OrderCalculable {
private:
    friend RegularSolid<Dodecahedron>;

    constexpr static double gold = (1 + std::sqrt(5.)) / 2;

    static void calculateStatic(const std::string &attr);

public:
    explicit Dodecahedron(const Matrix<3, 3> &orientation) : RegularSolid<Dodecahedron>{orientation} {};

    double projectionHalfsize(const Vector<3> &axis) const override;

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;
};


#endif //RSA3D_DODECAHEDRON_H
