//
// Created by PKua on 12.07.18.
//

#ifndef RSA3D_RHOMBICUBOCTAHEDRON_H
#define RSA3D_RHOMBICUBOCTAHEDRON_H


#include "RegularSolid.h"
#include "../../OrderCalculable.h"

class Rhombicuboctahedron : public RegularSolid<Rhombicuboctahedron>, public OrderCalculable {
private:
    friend RegularSolid<Rhombicuboctahedron>;

    static void calculateStatic(const std::string &attr);

public:
    explicit Rhombicuboctahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const override;

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;
};


#endif //RSA3D_RHOMBICUBOCTAHEDRON_H
