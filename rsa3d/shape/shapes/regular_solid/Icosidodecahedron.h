//
// Created by PKua on 13.07.18.
//

#ifndef RSA3D_ICOSIDODECAHEDRON_H
#define RSA3D_ICOSIDODECAHEDRON_H


#include "RegularSolid.h"

class Icosidodecahedron : public RegularSolid<Icosidodecahedron> {
private:
    friend RegularSolid<Icosidodecahedron>;

    constexpr static double goldRatio = (1 + std::sqrt(5.)) / 2;

    static void calculateStatic(const std::string &attr);

public:
    explicit Icosidodecahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const;         /* CRTP implement */
};


#endif //RSA3D_ICOSIDODECAHEDRON_H
