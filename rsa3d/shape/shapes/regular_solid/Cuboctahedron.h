//
// Created by PKua on 12.07.18.
//

#ifndef RSA3D_CUBEOCTAHEDRON_H
#define RSA3D_CUBEOCTAHEDRON_H


#include "RegularSolid.h"

class Cuboctahedron : public RegularSolid<Cuboctahedron> {
private:
    friend RegularSolid<Cuboctahedron>;

    static void calculateStatic(const std::string &attr);

public:
    explicit Cuboctahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const;         /* CRTP implement */
};


#endif //RSA3D_CUBEOCTAHEDRON_H
