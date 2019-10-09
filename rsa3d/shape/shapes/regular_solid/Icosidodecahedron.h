//
// Created by PKua on 13.07.18.
//

#ifndef RSA3D_ICOSIDODECAHEDRON_H
#define RSA3D_ICOSIDODECAHEDRON_H


#include "RegularSolid.h"
#include "Icosahedron.h"

class Icosidodecahedron : public RegularSolid<Icosidodecahedron> {
private:
    friend RegularSolid<Icosidodecahedron>;

    using SymmetryPlatonicSolid = Icosahedron;

    static ShapeData calculateStatic(const std::string &attr);

public:
    explicit Icosidodecahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const override;
};


#endif //RSA3D_ICOSIDODECAHEDRON_H
