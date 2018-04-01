//
// Created by PKua on 07.02.18.
//

#include "AnisotropicShape2D.h"

int AnisotropicShape2D::pointInside(BoundaryConditions *bc, double* position, double *orientation, double orientationRange){
	return this->pointInside(bc, position, orientation[0], orientation[0]+orientationRange);
}

Matrix<2, 2> AnisotropicShape2D::getRotationMatrix() const {
    return Matrix<2, 2>::rotation(this->getAngle());
}

Matrix<2, 2> AnisotropicShape2D::getAntiRotationMatrix() const {
    return Matrix<2, 2>::rotation(-this->getAngle());
}

// Keep angleFrom in [this->angle; this->angle + interval] range
//---------------------------------------------------------------------------------------------
void AnisotropicShape2D::normalizeAngleRange(double *angleFrom, double *angleTo, double interval) const {
    while (*angleFrom < this->getAngle()) {
        *angleFrom += interval;
        *angleTo += interval;
    }
    while (*angleFrom > this->getAngle() + interval) {
        *angleFrom -= interval;
        *angleTo -= interval;
    }
}

// Fix angle into [0; interval] range and return result
//---------------------------------------------------------------------------------------------
double AnisotropicShape2D::normalizeAngle(double angle, double interval) const {
    while (angle < 0)
        angle += interval;
    while (angle > interval)
        angle -= interval;
    return angle;
}

// This method is final and delegates to setAngle(double), so all deriving classes should
// override setAngle(double) method instead
//---------------------------------------------------------------------------------------------
void AnisotropicShape2D::rotate(double *v) {
    this->setAngle(this->getAngle() + *v);
}

// This method is final and delegates to setAngle(double), so all deriving classes should
// override setAngle(double) method instead
//---------------------------------------------------------------------------------------------
void AnisotropicShape2D::setOrientation(const double *orientation) {
    this->setAngle(*orientation);
}

// Keep angle in [0; getVoxelAngularSize()] range when setting.
//
// Derived classes should override this method when they want to keep track of the value of
// an angle, for example when storing exact vertices positions:
//
// void Derived::setAngle(double angle) {
//     AnisotropicShape2D::setAngle(angle);
//     this->calculateVertices();
// }
//---------------------------------------------------------------------------------------------
void AnisotropicShape2D::setAngle(double angle) {
    double interval = getVoxelAngularSize();
    double orientation = normalizeAngle(angle, interval);
    Shape::setOrientation(&orientation);  // Now use the original setter from Shape
}

double AnisotropicShape2D::getAngle() const{
    return this->getOrientation()[0];
}
