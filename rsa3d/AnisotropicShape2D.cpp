//
// Created by PKua on 07.02.18.
//

#include "AnisotropicShape2D.h"

double AnisotropicShape2D::getAngle() const{
	return this->orientation[0];
}

int AnisotropicShape2D::pointInside(BoundaryConditions *bc, double* position, double *orientation, double orientationRange){
	return this->pointInside(bc, position, orientation[0], orientation[0]+orientationRange);
}

Matrix<2, 2> AnisotropicShape2D::getRotationMatrix() const {
    return Matrix<2, 2>::rotation(this->orientation[0]);
}

Matrix<2, 2> AnisotropicShape2D::getAntiRotationMatrix() const {
    return Matrix<2, 2>::rotation(-this->orientation[0]);
}

// Keep angleFrom in [this->angle; this->angle + interval] range
//---------------------------------------------------------------------------------------------
void AnisotropicShape2D::normalizeAngleRange(double *angleFrom, double *angleTo, double interval) const {
    while (*angleFrom < this->orientation[0]) {
        *angleFrom += interval;
        *angleTo += interval;
    }
    while (*angleFrom > this->orientation[0] + interval) {
        *angleFrom -= interval;
        *angleTo -= interval;
    }
}
