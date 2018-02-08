//
// Created by PKua on 07.02.18.
//

#include <stdexcept>
#include "SpheroCylinder2D.h"

double SpheroCylinder2D::voxelSize;
double SpheroCylinder2D::neighbourListCellSize;
double SpheroCylinder2D::radius;
double SpheroCylinder2D::centerDistance;


void SpheroCylinder2D::initClass(const std::string &attr) {
    double ratio = std::stod(attr);
    if (ratio < 1) throw std::runtime_error("SpheroCylinder2D::initClass: ratio < 1");

    double normalization = std::sqrt(4 * (ratio - 1) + M_PI);
    radius = 1 / normalization;
    centerDistance = 2 * (ratio - 1) / normalization;
    neighbourListCellSize = centerDistance + 2 * radius;
    voxelSize = radius / M_SQRT2;
}

Shape<2> * SpheroCylinder2D::create(RND * rnd) {
    return (Shape<2>*)(new SpheroCylinder2D(rnd->nextValue() * M_PI));
}

double SpheroCylinder2D::getNeighbourListCellSize() {
    return neighbourListCellSize;
}

double SpheroCylinder2D::getVoxelSize() {
    return voxelSize;
}

double SpheroCylinder2D::getVolume() {
    return 1;
}

int SpheroCylinder2D::overlap(BoundaryConditions *bc, Shape *s) {
    SpheroCylinder2D other(*((SpheroCylinder2D*)s));
    double translation[2];
    bc->getTranslation(translation, this->getPosition(), s->getPosition());
    other.translate(translation);

    Vector<2> thisPos = this->getVectorPosition();
    Vector<2> otherPos = other.getVectorPosition();

    // Check bounding circle
    if (this->pointDistance2(otherPos) >= pow(2 * radius + centerDistance, 2))
        return false;

    // Check pallerogram inside
    Matrix<2, 2, double> thisRot = Matrix<2, 2>::rotation(this->angle);
    Matrix<2, 2, double> otherRot = Matrix<2, 2>::rotation(other.angle);
    Vector<2> otherPosAligned = thisRot.transpose() * (thisPos - otherPos);    // align x-axis to this capsule
    Vector<2> thisPosAligned = otherRot.transpose() * (thisPos - otherPos);    // align x-axis to colliding capsule
    if (fabs(otherPosAligned[1]) < fabs(centerDistance / 2 * sin(other.angle - this->angle)) &&
        fabs(thisPosAligned[1]) < fabs(centerDistance / 2 * sin(other.angle - this->angle)))
        return true;

    // Check 4 sides
    if (this->pointDistance2(thisPos + otherRot * Vector<2>{{centerDistance / 2 , 0}}, this->angle, otherPos) < 4 * radius * radius ||
        this->pointDistance2(thisPos - otherRot * Vector<2>{{centerDistance / 2 , 0}}, this->angle, otherPos) < 4 * radius * radius ||
        this->pointDistance2(thisPos + thisRot * Vector<2>{{centerDistance / 2 , 0}}, other.angle, otherPos) < 4 * radius * radius ||
        this->pointDistance2(thisPos - thisRot * Vector<2>{{centerDistance / 2 , 0}}, other.angle, otherPos) < 4 * radius * radius)
        return false;
}

double SpheroCylinder2D::pointDistance2(const Vector<2> &p) const
{
    return pointDistance2(this->getVectorPosition(), this->angle, p);
}

double SpheroCylinder2D::pointDistance2(const Vector<2> &arbitraryPos, double arbitraryAngle, const Vector<2> &p) const {
    Vector<2> diff = Matrix<2, 2>::rotation(-arbitraryAngle) * (arbitraryPos - p);

    if (fabs(diff[0]) < centerDistance / 2)    // middle (square) part
        return diff[1] * diff[1];
    else if (diff[0] < 0)     // left semicircle
        return (diff - Vector<2>{{-centerDistance / 2, 0}}).norm2();
    else                    // right semicircle
        return (diff - Vector<2>{{centerDistance / 2, 0}}).norm2();
}

int SpheroCylinder2D::pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) {
    return 0;
}

int SpheroCylinder2D::pointInside(BoundaryConditions *bc, double *da) {
    return AnisotropicShape2D::pointInside(bc, da);
}

std::string SpheroCylinder2D::toString() {
    return Shape::toString();
}

std::string SpheroCylinder2D::toPovray() const {
    return Shape::toPovray();
}

void SpheroCylinder2D::store(std::ostream &f) const {
    Shape::store(f);
}

void SpheroCylinder2D::restore(std::istream &f) {
    Shape::restore(f);
}

SpheroCylinder2D::SpheroCylinder2D(double angle) : angle(angle) {
}
