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
	void normalizeAngleRange(double *angleFrom, double *angleTo, double interval) const;
	double normalizeAngle(double angle, double interval) const;

public:
    // setOrientation(const double*) is delegated to setAngle(double)
    // rotate(double*) is delegated to setAngle(double)
    // pointInside(double* orientation, double orientationRange) -> pointInside(double angleFrom, double angleTo)
	void setOrientation(const double *orientation) final;
	void rotate(double *v) final;
    int pointInside(BoundaryConditions *bc, double* position, double *orientation, double orientationRange) final;

	double getAngle() const;
	virtual void setAngle(double angle);

    virtual int pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) = 0;
};


#endif //RSA3D_ANISOTROPICSHAPE2D_H
