//
// Created by PKua on 07.02.18.
//

#ifndef RSA3D_ANISOTROPICSHAPE2D_H
#define RSA3D_ANISOTROPICSHAPE2D_H


#include "Shape.h"
#include "Matrix.h"

class AnisotropicShape2D : public Shape<2, 1> {

protected:

    Matrix<2, 2> getAntiRotationMatrix() const;
    Matrix<2, 2> getRotationMatrix() const;

public:

    double getAngle() const;

    virtual int pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) = 0;
	virtual int pointInside(BoundaryConditions *bc, double* position, double *orientation, double orientationRange);

    void normalizeAngleRange(double angleFrom, double angleTo, double interval) const;
};


#endif //RSA3D_ANISOTROPICSHAPE2D_H
