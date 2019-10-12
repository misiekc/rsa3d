//
// Created by PKua on 14.07.18.
//

#ifndef RSA3D_TRUNCATEDICOSAHEDRON_H
#define RSA3D_TRUNCATEDICOSAHEDRON_H


#include "RegularSolid.h"
#include "Icosahedron.h"

class TruncatedIcosahedron : public RegularSolid<TruncatedIcosahedron> {
private:
    friend RegularSolid<TruncatedIcosahedron>;

    using SymmetryPlatonicSolid = Icosahedron;

    static ShapeData calculateStatic(const std::string &attr);

public:
    explicit TruncatedIcosahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const override;
};


#endif //RSA3D_TRUNCATEDICOSAHEDRON_H
