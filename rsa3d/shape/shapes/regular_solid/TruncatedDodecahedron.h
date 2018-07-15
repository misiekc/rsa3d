//
// Created by PKua on 13.07.18.
//

#ifndef RSA3D_TRUNCATEDDODECAHEDRON_H
#define RSA3D_TRUNCATEDDODECAHEDRON_H


#include "RegularSolid.h"

class TruncatedDodecahedron : public RegularSolid<TruncatedDodecahedron> {
private:
    friend RegularSolid<TruncatedDodecahedron>;

    constexpr static double goldRatio = (1 + std::sqrt(5.)) / 2;

    static void calculateStatic(const std::string &attr);

public:
    explicit TruncatedDodecahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const;         /* CRTP implement */
};


#endif //RSA3D_TRUNCATEDDODECAHEDRON_H