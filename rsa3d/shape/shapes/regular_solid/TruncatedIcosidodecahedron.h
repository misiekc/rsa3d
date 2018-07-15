//
// Created by PKua on 15.07.18.
//

#ifndef RSA3D_TRUNCATEDICOSIDODECAHEDRON_H
#define RSA3D_TRUNCATEDICOSIDODECAHEDRON_H


#include "RegularSolid.h"

class TruncatedIcosidodecahedron : public RegularSolid<TruncatedIcosidodecahedron> {
private:
    friend RegularSolid<TruncatedIcosidodecahedron>;

    constexpr static double goldRatio = (1 + std::sqrt(5.)) / 2;

    static void calculateStatic(const std::string &attr);

public:
    explicit TruncatedIcosidodecahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const;         /* CRTP implement */
};


#endif //RSA3D_TRUNCATEDICOSIDODECAHEDRON_H
