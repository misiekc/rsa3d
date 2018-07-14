//
// Created by PKua on 14.07.18.
//

#ifndef RSA3D_RHOMBICOSIDODECAHEDRON_H
#define RSA3D_RHOMBICOSIDODECAHEDRON_H


#include "RegularSolid.h"

class Rhombicosidodecahedron : public RegularSolid<Rhombicosidodecahedron> {
private:
    friend RegularSolid<Rhombicosidodecahedron>;

    constexpr static double goldRatio = (1 + std::sqrt(5.)) / 2;

    static void calculateStatic(const std::string &attr);

public:
    explicit Rhombicosidodecahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const;         /* CRTP implement */
};


#endif //RSA3D_RHOMBICOSIDODECAHEDRON_H
