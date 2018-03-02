//
// Created by PKua on 07.02.18.
//

#include "AnisotropicShape2D.h"

/**
 * Lazy pointInside - it uses angle-range pointInside with full angle
 * @param bc BoundaryConditions used
 * @param da point to check
 * @return 0 if not inside
 */
int AnisotropicShape2D::pointInside(BoundaryConditions *bc, double *da) {
    return pointInside(bc, da, 0, 2 * M_PI);
}

double AnisotropicShape2D::getAngle() const {
    return this->orientation[0];
}

void AnisotropicShape2D::setAngle(double angle) {
    this->orientation[0] = angle;
}

std::string AnisotropicShape2D::toWolfram() const {
    return "";
}

Matrix<2, 2> AnisotropicShape2D::getRotationMatrix() const {
    return Matrix<2, 2>::rotation(this->getAngle());
}

Matrix<2, 2> AnisotropicShape2D::getAntiRotationMatrix() const {
    return Matrix<2, 2>::rotation(-this->getAngle());
}

// Keep angleFrom in [this->angle; this->angle + interval] range
//---------------------------------------------------------------------------------------------
void AnisotropicShape2D::normalizeAngleRange(double &angleFrom, double &angleTo, double interval) const {
    while (angleFrom < this->getAngle()) {
        angleFrom += interval;
        angleTo += interval;
    }
    while (angleFrom > this->getAngle() + interval) {
        angleFrom -= interval;
        angleTo -= interval;
    }
}

void AnisotropicShape2D::rotate(double *v) {
    this->setAngle(this->getAngle() + v[0]);
}
