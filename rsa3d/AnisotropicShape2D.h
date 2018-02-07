//
// Created by PKua on 07.02.18.
//

#ifndef RSA3D_ANISOTROPICSHAPE2D_H
#define RSA3D_ANISOTROPICSHAPE2D_H


#include "Shape.h"

class AnisotropicShape2D : public Shape<2> {
public:
    virtual int pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) = 0;

    int pointInside(BoundaryConditions *bc, double *da) override;
};


#endif //RSA3D_ANISOTROPICSHAPE2D_H
