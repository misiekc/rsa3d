//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_TETRAHEDRON_H
#define RSA3D_TETRAHEDRON_H


#include "RegularSolid.h"

class Tetrahedron : public RegularSolid<Tetrahedron> {
private:
    friend RegularSolid<Tetrahedron>;

    static void calculateStatic(const std::string &attr);

public:
    explicit Tetrahedron(const Matrix<3, 3> &orientation) : RegularSolid<Tetrahedron>{orientation} {};

    double projectionHalfsize(const Vector<3> &axis) const;     /* CRTP implement */
    bool isSeparatingAxis(const Vector<3> &axis, const Tetrahedron &other,
                          const Vector<3> &distance) const;     /* CRTP override */

    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     double orientationRange) const override;

};


#endif //RSA3D_TETRAHEDRON_H
