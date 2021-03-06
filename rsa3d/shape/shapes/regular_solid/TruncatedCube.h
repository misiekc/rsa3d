//
// Created by PKua on 12.07.18.
//

#ifndef RSA3D_TRUNCATEDCUBE_H
#define RSA3D_TRUNCATEDCUBE_H


#include "RegularSolid.h"
#include "Octahedron.h"

class TruncatedCube : public RegularSolid<TruncatedCube> {
private:
    friend RegularSolid<TruncatedCube>;

    using SymmetryPlatonicSolid = Octahedron;

    static ShapeData calculateStatic(const std::string &attr);

public:
    explicit TruncatedCube(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const override;
};


#endif //RSA3D_TRUNCATEDCUBE_H
