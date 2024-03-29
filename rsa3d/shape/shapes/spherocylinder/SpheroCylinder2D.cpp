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
ShapeStaticInfo<2, 1> SpheroCylinder2D::spherocylinderShapeInfo;

void SpheroCylinder2D::calculateStatic(const std::string &attr) {
    std::istringstream attrStream(attr);
    double firstArg, secondArg;
    attrStream >> firstArg;
    ValidateMsg(attrStream, "Invalid attibutes format. Expecting: [ratio] OR [distance] [radius]");
    attrStream >> secondArg;
    bool hasSecondArg = static_cast<bool>(attrStream);

    if (hasSecondArg) {
        // arguments: [distance between disks] [radius]

        halfDistance = firstArg/2;
        radius = secondArg;
        Validate(halfDistance >= 0);
        Validate(radius > 0);
    } else {
        // arguments: [ratio]; volume normalize to 1

        double ratio = firstArg;
        Validate(ratio >= 1);

        double normalization = std::sqrt(4 * (ratio - 1) + M_PI);
        radius = 1 / normalization;
        halfDistance = (ratio - 1) / normalization;
    }

    centerVector = Vector<2>{{halfDistance , 0}};
}

void SpheroCylinder2D::initClass(const std::string &attr) {
    SpheroCylinder2D::calculateStatic(attr);

    // First of all, initialize spherocylinder own ShapeInfo.
    // It won't be lost after diskopolygon-type shapes rerun ShapeFactory::initShapeClass (because is used
    // SpheroCylinder2D internally, it has to initialize it at the beginning
    SpheroCylinder2D::spherocylinderShapeInfo.setCircumsphereRadius(halfDistance + radius);
    SpheroCylinder2D::spherocylinderShapeInfo.setInsphereRadius(radius);
    SpheroCylinder2D::spherocylinderShapeInfo.setAngularVoxelSize(M_PI);
    SpheroCylinder2D::spherocylinderShapeInfo.setSupportsSaturation(true);
    SpheroCylinder2D::spherocylinderShapeInfo.setDefaultCreateShapeImpl <SpheroCylinder2D> ();

    // Set normal, global shape info to that one.
    // This one will be overwritten in case of diskopolygon simulation, but in that case only
    // SpheroCylinder2D::overlap and SpheroCylinder2D::pointInside methods matter, and they are not affected by the
    // global Shape's static info
    Shape::setShapeStaticInfo(SpheroCylinder2D::spherocylinderShapeInfo);
}

double SpheroCylinder2D::getVolume(unsigned short dim) const {
	switch (dim){
	case 2:
		return 4*SpheroCylinder2D::halfDistance*SpheroCylinder2D::radius +
				M_PI*SpheroCylinder2D::radius*SpheroCylinder2D::radius;
		break;
	case 1:
	{
		double diameter = 2*SpheroCylinder2D::radius;
		double distance = 2*SpheroCylinder2D::halfDistance;
		double angle = this->getAngle();
		double vertexAngle = std::atan(diameter / distance);
	    if (angle < vertexAngle)
	        return distance * std::cos(angle) + std::sqrt(diameter * diameter - distance * distance * std::sin(angle) * std::sin(angle));
	    else if (angle < M_PI - vertexAngle)
	        return diameter / std::sin(angle);
	    else if (angle < M_PI + vertexAngle)
	        return -distance * std::cos(angle) + std::sqrt(diameter * diameter - distance * distance * std::sin(angle) * std::sin(angle));
	    else if (angle < 2 * M_PI - vertexAngle)
	        return -diameter / std::sin(angle);
//	    else
	        return distance * std::cos(angle) + std::sqrt(diameter * diameter - distance * distance * std::sin(angle) * std::sin(angle));
		break;
	}
	default:
        throw std::runtime_error ("SpheroCylinder2D supports only 2D and 1D packings");
		break;
	}
}

bool SpheroCylinder2D::overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const {
    // Use "local", SpheroCylinder2D ShapeStaticInfo, because for diskopolygon classes, the global one is overriden
    switch (this->overlapEarlyRejection(bc, s, SpheroCylinder2D::spherocylinderShapeInfo)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }

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
    AnisotropicShape2D::normalizeAngleRange(this->getAngle(), &angleFrom, &angleTo, M_PI);

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
	out << "0.0>" << std::endl << "	texture { finish { ambient 1 diffuse 0 } pigment { color Red} }" << std::endl << "}" << std::endl;

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

void SpheroCylinder2D::setAngle(double angle) {
    angle = this->normalizeAngle(angle, M_PI);
    AnisotropicShape2D::setAngle(angle);
}
