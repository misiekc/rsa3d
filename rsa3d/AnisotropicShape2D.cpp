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

AnisotropicShape2D::AnisotropicShape2D(double angle) : angle(angle) {

}

double AnisotropicShape2D::getAngle() const {
    return angle;
}

void AnisotropicShape2D::setAngle(double angle) {
    AnisotropicShape2D::angle = angle;
}

std::string AnisotropicShape2D::toWolfram() const {
    return "";
}

Matrix<2, 2> AnisotropicShape2D::getRotationMatrix() const {
    return Matrix<2, 2>::rotation(angle);
}

Matrix<2, 2> AnisotropicShape2D::getAntiRotationMatrix() const {
    return Matrix<2, 2>::rotation(-angle);
}