//
// Created by PKua on 07.02.18.
//

#include "AnisotropicShape2D.h"

bool AnisotropicShape2D::pointInside(BoundaryConditions *bc, double* position, const std::array<double, 1> &orientation,
                                    double orientationRange) const {
	return this->pointInside(bc, position, orientation[0], orientation[0]+orientationRange);
}

bool AnisotropicShape2D::pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) const{
	return 0;
}

Matrix<2, 2> AnisotropicShape2D::getRotationMatrix() const {
    return Matrix<2, 2>::rotation(this->getAngle());
}

Matrix<2, 2> AnisotropicShape2D::getAntiRotationMatrix() const {
    return Matrix<2, 2>::rotation(-this->getAngle());
}

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

double AnisotropicShape2D::normalizeAngle(double angle, double interval) const {
    while (angle < 0)
        angle += interval;
    while (angle > interval)
        angle -= interval;
    return angle;
}

void AnisotropicShape2D::setOrientation(const double *orientation) {
    this->setAngle(*orientation);
}

void AnisotropicShape2D::setAngle(double angle) {
    double interval = getVoxelAngularSize();
    double orientation = normalizeAngle(angle, interval);
    Shape::setOrientation(&orientation);  // Now use the original setter from Shape
}

double AnisotropicShape2D::getAngle() const{
    return this->getOrientation()[0];
}
