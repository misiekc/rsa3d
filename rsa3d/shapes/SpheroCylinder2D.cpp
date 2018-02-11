//
// Created by PKua on 07.02.18.
//

#include <stdexcept>
#include <iomanip>
#include "SpheroCylinder2D.h"

#define EPSILON 0.0000001

double SpheroCylinder2D::voxelSize;
double SpheroCylinder2D::neighbourListCellSize;
double SpheroCylinder2D::radius;
double SpheroCylinder2D::halfDistance;
Vector<2> SpheroCylinder2D::centerVector;


void SpheroCylinder2D::initClass(const std::string &attr) {
    double ratio = std::stod(attr);
    if (ratio < 1) throw std::runtime_error("SpheroCylinder2D::initClass: ratio < 1");

    double normalization = std::sqrt(4 * (ratio - 1) + M_PI);
    radius = 1 / normalization;
    halfDistance = (ratio - 1) / normalization;
    neighbourListCellSize = (halfDistance + radius) * 2;
    voxelSize = radius / M_SQRT2;
    centerVector = Vector<2>{{halfDistance , 0}};
}

Shape<2> * SpheroCylinder2D::create(RND * rnd) {
    return new SpheroCylinder2D(rnd->nextValue() * M_PI);
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
    this->applyBC(bc, &other);
    return withinExclusionZone(other.getVectorPosition(), other.angle);
}

Matrix<2, 2> SpheroCylinder2D::getRotationMatrix() const {
    return Matrix<2, 2>::rotation(angle);
}

Matrix<2, 2> SpheroCylinder2D::getAntiRotationMatrix() const {
    return Matrix<2, 2>::rotation(-angle);
}

double SpheroCylinder2D::pointDistance2(const Vector<2> &p) const
{
    return pointDistance2(this->getVectorPosition(), this->angle, p);
}

double SpheroCylinder2D::pointDistance2(const Vector<2> &pos, double angle, const Vector<2> &point) {
    Vector<2> diff = Matrix<2, 2>::rotation(-angle) * (pos - point);

    if (fabs(diff[0]) < halfDistance)    // middle (square) part
        return diff[1] * diff[1];
    else if (diff[0] < 0)   // left semicircle
        return (diff + centerVector).norm2();
    else                    // right semicircle
        return (diff - centerVector).norm2();
}

int SpheroCylinder2D::pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) {
    this->normalizeRange(angleFrom, angleTo);

    // If angleTo not in normal range (see normalizeRange), divide it and check separately
    if (angleTo > this->angle + M_PI) {
        return pointInside(bc, da, angleFrom, this->angle + M_PI - EPSILON) &&
               pointInside(bc, da, this->angle, angleTo - M_PI);
    }

    Vector<2> pointPos = applyBC(bc, da);
    Vector<2> pointPosThisAligned = getAntiRotationMatrix() * (pointPos - this->getVectorPosition());
    double angleFromAligned = angleFrom - this->angle;
    double angleToAligned = angleTo - this->angle;
    double pointAngleToOrigin = getAngleToOrigin(pointPosThisAligned);
    double pointAngleToRDisk = getAngleToOrigin(pointPosThisAligned - centerVector);
    double pointAngleToLDisk = getAngleToOrigin(pointPosThisAligned + centerVector);

    // Check "special areas" around disks
    if (angleInRange(pointAngleToOrigin, 0, M_PI/2)) {      // first quarter (right disk)
        if (angleInRange(pointAngleToRDisk, angleFromAligned - M_PI/2, angleToAligned - M_PI/2))
            return this->pointDistance2(pointPos) < 4 * radius * radius;
    } else if (angleInRange(pointAngleToOrigin, M_PI/2, 3*M_PI/2)) {    // second and third (left disk)
        if (angleInRange(pointAngleToLDisk, angleFromAligned + M_PI/2, angleToAligned + M_PI/2))
            return this->pointDistance2(pointPos) < 4 * radius * radius;
    } else if (angleInRange(pointAngleToOrigin, 3*M_PI/2, 2*M_PI)) {    // fourth (right disk)
        if (angleInRange(pointAngleToRDisk, angleFromAligned + 3*M_PI/2, angleToAligned + 3*M_PI/2))
            return this->pointDistance2(pointPos) < 4 * radius * radius;
    }

    // Check "standard area"
    return this->withinExclusionZone(pointPos, angleFrom) && this->withinExclusionZone(pointPos, angleTo);
}

// Keep angleFrom in [this->angle; this->angle + pi] range
//---------------------------------------------------------------------------------------------
void SpheroCylinder2D::normalizeRange(double &angleFrom, double &angleTo) const {
    while (angleFrom < angle) {
        angleFrom += M_PI;
        angleTo += M_PI;
    }
    while (angleFrom > angle + M_PI) {
        angleFrom -= M_PI;
        angleTo -= M_PI;
    }
}

bool SpheroCylinder2D::angleInRange(double angle, double rangeStart, double rangeEnd) const {
    return rangeStart < angle && angle <= rangeEnd;
}

int SpheroCylinder2D::pointInside(BoundaryConditions *bc, double *da) {
    double translationArr[2];
    Vector<2> translationVector(bc->getTranslation(translationArr, getPosition(), da));
    Vector<2> point = Vector<2>(da) + translationVector;
    return pointDistance2(point) < 4 * radius * radius;
}

std::string SpheroCylinder2D::toString() {
    std::stringstream out;
    out << "SpheroCylinder2D{radius: " << radius << "; halfDistance: " << halfDistance;
    out << "; pos: " << this->getVectorPosition() << "; angle: " << angle << "}";
    return out.str();
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

std::string SpheroCylinder2D::toWolfram() const {
    std::stringstream out;

    out << std::fixed;
    out << "GeometricTransformation[{Rectangle[{-" << halfDistance << ", -" << radius << "}, {" << halfDistance << ", " << radius << "}]," << std::endl;
    out << "    Disk[{-" << halfDistance << ", 0}, " << radius << "]," << std::endl;
    out << "    Disk[{" << halfDistance << ", 0}, " << radius << "]}," << std::endl;
    out << "    {RotationMatrix[" << this->angle << "], " << this->getVectorPosition() << "}]";

    return out.str();
}

double SpheroCylinder2D::getAngleToOrigin(const Vector<2> & point) {
    double angle = atan2(point[1], point[0]);
    if(angle < 0)
        return angle + 2 * M_PI;
    else
        return angle;
}

bool SpheroCylinder2D::withinExclusionZone(const Vector<2> &pointPos, double angle) {
    Vector<2> thisPos = this->getVectorPosition();

    // Check bounding spherocylinder
    if (this->pointDistance2(pointPos) >= pow(2 * (radius + halfDistance), 2))
        return false;

    // Check small parallelogram inside
    Matrix<2, 2> thisRot = this->getRotationMatrix();
    Matrix<2, 2> pointRot = Matrix<2, 2>::rotation(angle);
    Vector<2> pointPosThisAligned = thisRot.transpose() * (pointPos - thisPos);     // align x-axis to this capsule
    Vector<2> thisPosPointAligned = pointRot.transpose() * (thisPos - pointPos);    // align x-axis to colliding capsule
    if (fabs(pointPosThisAligned[1]) < fabs(halfDistance * sin(angle - this->angle)) &&
        fabs(thisPosPointAligned[1]) < fabs(halfDistance * sin(this->angle - angle)))
        return true;

    // Check big parallelogram with round vertices
    return this->pointDistance2(thisPos + pointRot * centerVector, this->angle, pointPos) < 4 * radius * radius ||
           this->pointDistance2(thisPos - pointRot * centerVector, this->angle, pointPos) < 4 * radius * radius ||
           this->pointDistance2(thisPos + thisRot * centerVector, angle, pointPos) < 4 * radius * radius ||
           this->pointDistance2(thisPos - thisRot * centerVector, angle, pointPos) < 4 * radius * radius;
}
