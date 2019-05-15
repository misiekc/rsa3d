//
// Created by PKua on 07.02.18.
//

#include <stdexcept>
#include <iomanip>
#include "SpheroCylinder2D.h"
#include "../../../utils/Assertions.h"

static const double EPSILON = 0.0000000001;

double SpheroCylinder2D::radius;
double SpheroCylinder2D::halfDistance;
Vector<2> SpheroCylinder2D::centerVector;

void SpheroCylinder2D::calculateStatic(const std::string &attr) {
    double ratio = std::stod(attr);
    ValidateMsg(ratio >= 1, "Ratio has to be >= 1");

    double normalization = std::sqrt(4 * (ratio - 1) + M_PI);
    radius = 1 / normalization;
    halfDistance = (ratio - 1) / normalization;
    centerVector = Vector<2>{{halfDistance , 0}};
}

void SpheroCylinder2D::initClass(const std::string &attr) {
    SpheroCylinder2D::calculateStatic(attr);

    Shape::setNeighbourListCellSize((halfDistance + radius) * 2);
    Shape::setVoxelSpatialSize(M_SQRT2 * radius);
    Shape::setVoxelAngularSize(M_PI);
	Shape::setSupportsSaturation(true);
    Shape::setDefaultCreateShapeImpl <SpheroCylinder2D> ();
}

double SpheroCylinder2D::getVolume() const {
    return 1;
}

bool SpheroCylinder2D::overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const {
    SpheroCylinder2D other = dynamic_cast<const SpheroCylinder2D&>(*s);
    this->applyBC(bc, &other);
    return withinExclusionZone(other.getPosition(), other.getAngle());
}

double SpheroCylinder2D::pointDistance2(const Vector<2> &p) const
{
    return pointDistance2(this->getPosition(), this->getAngle(), p);
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

bool SpheroCylinder2D::pointInside(BoundaryConditions<2> *bc, const Vector<2> &da, double angleFrom,
                                   double angleTo) const {
    this->normalizeAngleRange(&angleFrom, &angleTo, M_PI);

    // If angleTo not in normal range (see normalizeAngleRange), divide it and check separately
    if (angleTo > this->getAngle() + M_PI) {
        return pointInside(bc, da, angleFrom, this->getAngle() + M_PI - EPSILON) &&
               pointInside(bc, da, this->getAngle(), angleTo - M_PI);
    }

    Vector<2> translationVector = bc->getTranslation(this->getPosition(), da);
    Vector<2> thisPos(this->getPosition());
    Vector<2> pointPos = da + translationVector;
    Vector<2> pointPosThisAligned = getAntiRotationMatrix() * (pointPos - thisPos);
    double angleFromAligned = angleFrom - this->getAngle();
    double angleToAligned = angleTo - this->getAngle();
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

bool SpheroCylinder2D::angleInRange(double angle, double rangeStart, double rangeEnd) const {
    return rangeStart < angle && angle <= rangeEnd;
}

/*
int SpheroCylinder2D::pointInside(BoundaryConditions *bc, double *da) {
    double translationArr[2];
    Vector<2> translationVector(bc->getTranslation(translationArr, getPosition(), da));
    Vector<2> point = Vector<2>(da) + translationVector;
    return pointDistance2(point) < 4 * radius * radius;
}
*/

std::string SpheroCylinder2D::toString() const {
    std::stringstream out;
    Vector<2> thisPos(this->getPosition());
    out << "SpheroCylinder2D{radius: " << radius << "; halfDistance: " << halfDistance;
    out << "; pos: " << thisPos << "; angle: " << this->getAngle() * 180 / M_PI << "}";
    return out.str();
}

std::string SpheroCylinder2D::toPovray() const {
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);
	out << "  union {" << std::endl;;

	out << "    disc { < " << -halfDistance << ", 0.0, 0.0002>, <0.0, 0.0, 1.0>, " << radius << " }" << std::endl;
	out << "    disc { < " <<  halfDistance << ", 0.0, 0.0002>, <0.0, 0.0, 1.0>, " << radius << " }" << std::endl;
	out << "    polygon { 5, < " << -halfDistance << ", " << -radius << ", 0.0002>, < " << -halfDistance << ", " << radius << ", 0.0002>, < " << halfDistance << ", " << radius << ", 0.0002>, < " << halfDistance << ", " << -radius << ", 0.0002> < " << -halfDistance << ", " << -radius << ", 0.0002> }" << std::endl;
	out << "    rotate <0, 0, " << (180*this->getAngle()/M_PI) << ">" << std::endl;
	out << "	translate <";

    Vector<2> position = this->getPosition();
	for(unsigned short i=0; i<2; i++)
		out << (position[i]) << ", ";
	out << "0.0>" << std::endl << "	texture { pigment { color Red } }" << std::endl << "}" << std::endl;

	return out.str();

}

void SpheroCylinder2D::store(std::ostream &f) const {
    Shape::store(f);
}

void SpheroCylinder2D::restore(std::istream &f) {
    Shape::restore(f);
}

std::string SpheroCylinder2D::toWolfram() const {
    std::stringstream out;
    Vector<2> thisPos(this->getPosition());

    out << std::fixed;
    out << "GeometricTransformation[{Rectangle[{-" << halfDistance << ", -" << radius << "}, {" << halfDistance << ", " << radius << "}]," << std::endl;
    out << "    Disk[{-" << halfDistance << ", 0}, " << radius << "]," << std::endl;
    out << "    Disk[{" << halfDistance << ", 0}, " << radius << "]}," << std::endl;
    out << "    {RotationMatrix[" << this->getAngle() << "], " << thisPos << "}]";

    return out.str();
}

bool SpheroCylinder2D::withinExclusionZone(const Vector<2> &pointPos, double angle) const {
    Vector<2> thisPos(this->getPosition());

    // Check bounding spherocylinder
    if (this->pointDistance2(pointPos) >= pow(2 * (radius + halfDistance), 2))
        return false;

    // Check small parallelogram inside
    Matrix<2, 2> thisRot = this->getRotationMatrix();
    Matrix<2, 2> pointRot = Matrix<2, 2>::rotation(angle);
    Vector<2> pointPosThisAligned = thisRot.transpose() * (pointPos - thisPos);     // align x-axis to this capsule
    Vector<2> thisPosPointAligned = pointRot.transpose() * (thisPos - pointPos);    // align x-axis to colliding capsule
    if (fabs(pointPosThisAligned[1]) < fabs(halfDistance * sin(angle - this->getAngle() )) &&
        fabs(thisPosPointAligned[1]) < fabs(halfDistance * sin(this->getAngle() - angle)))
        return true;

    // Check big parallelogram with round vertices
    return this->pointDistance2(thisPos + pointRot * centerVector, this->getAngle(), pointPos) < 4 * radius * radius ||
           this->pointDistance2(thisPos - pointRot * centerVector, this->getAngle(), pointPos) < 4 * radius * radius ||
           this->pointDistance2(thisPos + thisRot * centerVector, angle, pointPos) < 4 * radius * radius ||
           this->pointDistance2(thisPos - thisRot * centerVector, angle, pointPos) < 4 * radius * radius;
}

Shape<2, 1> *SpheroCylinder2D::clone() const {
    return new SpheroCylinder2D(*this);
}
