//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_DODECAHEDRON_H
#define RSA3D_DODECAHEDRON_H


#include "RegularSolid.h"
#include "Icosahedron.h"

class Dodecahedron : public RegularSolid<Dodecahedron> {
private:
    friend RegularSolid<Dodecahedron>;

    using SymmetryPlatonicSolid = Icosahedron;

    static ShapeData calculateStatic(const std::string &attr);

public:
    explicit Dodecahedron(const Matrix<3, 3> &orientation) : RegularSolid<Dodecahedron>{orientation} {};

    double projectionHalfsize(const Vector<3> &axis) const override;
};


#endif //RSA3D_DODECAHEDRON_H
