//
// Created by PKua on 07.02.18.
//

#ifndef RSA3D_ANISOTROPICSHAPE2D_H
#define RSA3D_ANISOTROPICSHAPE2D_H


#include "Shape.h"
#include "AnisotropicShape2D.h"

class AnisotropicShape2D : public Shape<2> {

protected:
    double angle{};

    Matrix<2, 2> getAntiRotationMatrix() const;
    Matrix<2, 2> getRotationMatrix() const;

public:
    AnisotropicShape2D() = default;
    AnisotropicShape2D(double angle);

    virtual int pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) = 0;
    int pointInside(BoundaryConditions *bc, double *da) override;
    virtual double getAngle() const;
    virtual void setAngle(double angle);
    virtual std::string toWolfram() const;
};


#endif //RSA3D_ANISOTROPICSHAPE2D_H
