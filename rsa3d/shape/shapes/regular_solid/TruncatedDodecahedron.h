//
// Created by PKua on 13.07.18.
//

#ifndef RSA3D_TRUNCATEDDODECAHEDRON_H
#define RSA3D_TRUNCATEDDODECAHEDRON_H


#include "RegularSolid.h"
#include "../../OrderCalculable.h"

class TruncatedDodecahedron : public RegularSolid<TruncatedDodecahedron> {
private:
    friend RegularSolid<TruncatedDodecahedron>;

    static void calculateStatic(const std::string &attr);

public:
    explicit TruncatedDodecahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const override;
};


#endif //RSA3D_TRUNCATEDDODECAHEDRON_H
