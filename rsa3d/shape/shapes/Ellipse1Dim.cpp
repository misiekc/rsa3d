/*
 * Ellipse.cpp
 *
 *  Created on: 10.08.2017
 *      Author: ciesla
 */

#include "Ellipse1Dim.h"
#include "../../utils/Utils.h"
#include "../../geometry/Geometry.h"


static const double EPSILON = 0.0000000001;

double Ellipse1Dim::longSemiAxis;
double Ellipse1Dim::shortSemiAxis;

void Ellipse1Dim::calculateU() {
	this->u[0]  = cos(this->getAngle()); this->u[1]  = sin(this->getAngle());
	this->uT[0] = -sin(this->getAngle()); this->uT[1] =  cos(this->getAngle());
}

Ellipse1Dim::Ellipse1Dim() : ConvexShape<1,1>(){
	this->a = Ellipse1Dim::longSemiAxis;
	this->b = Ellipse1Dim::shortSemiAxis;
	this->setAngle(0);
}

void Ellipse1Dim::initClass(const std::string &args){
	double ratio = std::stod(args);
	Ellipse1Dim::shortSemiAxis = sqrt(1.0/(M_PI*ratio));
	Ellipse1Dim::longSemiAxis = ratio*shortSemiAxis;
	
	ShapeStaticInfo<1, 1> shapeInfo;
	shapeInfo.setCircumsphereRadius(longSemiAxis);
	shapeInfo.setInsphereRadius(shortSemiAxis);
	shapeInfo.setAngularVoxelSize({M_PI});
	shapeInfo.setSupportsSaturation(true);
	shapeInfo.setDefaultCreateShapeImpl <Ellipse1Dim> ();
	
	Shape::setShapeStaticInfo(shapeInfo);
}

double Ellipse1Dim::calculateF(double* r, double g) const {
	double d1 = (r[0]*this->u[0] + r[1]*this->u[1])  /this->a;
	double d2 = (r[0]*this->uT[0] + r[1]*this->uT[1]) /this->b;

	return 1 + g - d1*d1 - d2*d2;
}

void Ellipse1Dim::setAngle(double angle){
	Orientation<1> interval = getAngularVoxelSize();
	Orientation<1> orientation{{normalizeAngle(angle, interval[0])}};
	Shape::setOrientation(orientation);  // Now use the original setter from Shape

	this->calculateU();
}

bool Ellipse1Dim::overlap(BoundaryConditions<1> *bc, const Shape<1, 1> *s) const {
    switch (this->overlapEarlyRejection(bc, s)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }
    
	Ellipse1Dim es = dynamic_cast<const Ellipse1Dim &>(*s);
    this->applyBC(bc, &es);
	double d;
	d = (this->a/this->b - this->b/this->a)*sin(this->getAngle() - es.getAngle());
	double g = 2 + d*d;
	double r[] = {this->getPosition()[0] - es.getPosition()[0], 0};
	double f1 = this->calculateF(r, g);
	double f2 = es.calculateF(r, g);
	double psi = 4*(f1*f1-3*f2)*(f2*f2-3*f1) - (9-f1*f2)*(9-f1*f2);
	if (psi>0 && (f1<0 || f2<0)){
		return false;
	}
	return true;
}

double Ellipse1Dim::getVolume(unsigned short dim) const {
    if (dim != 1)
        throw std::runtime_error ("Ellipse1Dim supports only 1D packings");

    return M_PI*this->a*this->b;
}

// TODO not working properly

/*bool Ellipse1Dim::pointInside(BoundaryConditions *bc, double* da) const {
	double ta[2];
	double tmp[2];

	double *position = this->getPosition();

	tmp[0] = da[0]; tmp[1] = da[1];
	bc->getTranslation(ta, position, tmp);
	da[0] += ta[0];
	da[1] += ta[1];

	da[0] -= position[0];
	da[1] -= position[1];

	rotate2D(da, -this->getAngle());

	double dx = da[0]/(this->a+this->b);
	double dy = da[1]/(2*this->b);
	return (dx*dx+dy*dy < 1);
}*/

bool Ellipse1Dim::pointInside(BoundaryConditions<1> *bc, const Vector<1> &other, double angleFrom, double angleTo) const
{
    // Check exlusion zones for angle interval endpoinds - angleFrom and angleTo
    Ellipse1Dim ellTmp;
    ellTmp.translate(other);

    ellTmp.setAngle(angleFrom);
    if (!this->overlap(bc, &ellTmp))
        return false;
    ellTmp.setAngle(angleTo);
    if (!this->overlap(bc, &ellTmp))
        return false;

    // Now check the "special areas"
	return pointInsideSpecialArea(bc, other, angleFrom, angleTo);
}

bool Ellipse1Dim::pointInsideSpecialArea(BoundaryConditions<1> *bc, const Vector<1> &other, double angleFrom,
									double angleTo) const {

	this->normalizeAngleRange(&angleFrom, &angleTo, M_PI);
	// now angleFrom is in [this->getAngle(), this->getAngle() + M_PI)

    // If angleTo is not in not in [this->getAngle(), this->getAngle() + M_PI), then two checks are made separately
	if (angleTo > this->getAngle() + M_PI) {
        return pointInsideSpecialArea(bc, other, angleFrom, this->getAngle() + M_PI - EPSILON) &&
               pointInsideSpecialArea(bc, other, this->getAngle(), angleTo - M_PI);
    }

    // Now angleTo - angleFrom <= pi. Align ellipse with coordinate system. Transform angles' range
    // (and divide if necessary) so that its fully contained in [0, pi] range
    Vector<1> translationVector = bc->getTranslation(this->getPosition(), other);
    Vector<1> thisPos(this->getPosition());
    Vector<1> otherPos = Vector<1>(other) + translationVector;
    Vector<2> tmp({thisPos[0] - otherPos[0], 0});
	Vector<2> otherAligned = getAntiRotationMatrix() * (tmp);
	if (angleTo < this->getAngle())
		return pointInsideUnrotated(otherAligned, angleFrom - this->getAngle() + M_PI, angleTo - this->getAngle() + M_PI);
	else if (angleFrom < this->getAngle())
		return pointInsideUnrotated(otherAligned, angleFrom - this->getAngle() + M_PI, M_PI) &&
               pointInsideUnrotated(otherAligned, 0, angleTo - this->getAngle());
	else    // angleTo >= angle && angleFrom >= angle
		return pointInsideUnrotated(otherAligned, angleFrom - this->getAngle(), angleTo - this->getAngle());
}

bool Ellipse1Dim::pointInsideUnrotated(const Vector<2> &p, double angleFrom, double angleTo) const
{
	if (angleFrom <= M_PI/2 && angleTo >= M_PI/2)
    {
		if (p[0] <= 0 && p[1] >= 0 && withinAngle(p, angleFrom + M_PI/2, M_PI))
			return withinAngleCheckCollision(p, angleFrom + M_PI/2, M_PI);
		if (p[0] <= 0 && p[1] <= 0 && withinAngle(p, M_PI + EPSILON, angleTo + M_PI/2))
			return withinAngleCheckCollision(p, M_PI + EPSILON, angleTo + M_PI/2);
		if (p[0] >= 0 && p[1] <= 0 && withinAngle(-p, angleFrom + M_PI/2, M_PI))
			return withinAngleCheckCollision(-p, angleFrom + M_PI/2, M_PI);
		if (p[0] >= 0 && p[1] >= 0 && withinAngle(-p, M_PI + EPSILON, angleTo + M_PI/2))
			return withinAngleCheckCollision(-p, M_PI + EPSILON, angleTo + M_PI/2);
	} 
    else if (angleFrom <= M_PI/2 && angleTo <= M_PI/2)
    {
		if (p[0] <= 0 && p[1] >= 0 && withinAngle(p, angleFrom + M_PI/2, angleTo + M_PI/2))
			return withinAngleCheckCollision(p, angleFrom + M_PI/2, angleTo + M_PI/2);
		if (p[0] >= 0 && p[1] <= 0 && withinAngle(-p, angleFrom + M_PI/2, angleTo + M_PI/2))
			return withinAngleCheckCollision(-p, angleFrom + M_PI/2, angleTo + M_PI/2);
	} 
    else // angleFrom >= M_PI/2 && angleTo >= M_PI/2
    {
		if (p[0] <= 0 && p[1] <= 0 && withinAngle(p, angleFrom + M_PI/2, angleTo + M_PI/2))
			return withinAngleCheckCollision(p, angleFrom + M_PI/2, angleTo + M_PI/2);
		if (p[0] >= 0 && p[1] >= 0 && withinAngle(-p, angleFrom + M_PI/2, angleTo + M_PI/2))
			return withinAngleCheckCollision(-p, angleFrom + M_PI/2, angleTo + M_PI/2);
	}

	return true;
}

bool Ellipse1Dim::withinAngle(const Vector<2> &p, double angleFrom, double angleTo) const
{
	double minT = atan2(b * sin(angleFrom), a * cos(angleFrom));
	double maxT = atan2(b * sin(angleTo), a * cos(angleTo));

	Vector<2> edgeMin{{cos(minT) * a, sin(minT) * b}};
	Vector<2> edgeMax{{cos(maxT) * a, sin(maxT) * b}};
	Vector<2> edgeCross = intersection::line_line(edgeMin, angleFrom, edgeMax, angleTo);
	double zoneAngle = getAngleToOrigin(p - edgeCross);
    
    return zoneAngle >= angleFrom && zoneAngle <= angleTo;
}

bool Ellipse1Dim::withinAngleCheckCollision(const Vector<2> &p, double lowerAngle, double upperAngle) const
{
	double minT = atan2(b * sin(lowerAngle), a * cos(lowerAngle));
	double maxT = atan2(b * sin(upperAngle), a * cos(upperAngle));

	return circleCollision(p, minT, maxT);
}

bool Ellipse1Dim::circleCollision(const Vector<2> &p, double tMin, double tMax) const
{
    Vector<2> start{{cos(tMin) * a, sin(tMin) * b}};
    Vector<2> end{{cos(tMax) * a, sin(tMax) * b}};
    if (collision::line_circle(p, b, start, end))
        return true;

    double slope1 = atan2(a * sin(tMin), b * cos(tMin)) + M_PI/2;
    double slope2 = atan2(a * sin(tMax), b * cos(tMax)) + M_PI/2;
    Vector<2> intersect = intersection::line_line(start, slope1, end, slope2);
	if (!collision::line_circle(p, b, start, intersect) && !collision::line_circle(p, b, intersect, end))
		return false;

	return circleCollision(p, tMin, (tMin + tMax) / 2) || circleCollision(p, (tMin + tMax) / 2, tMax);
}

std::string Ellipse1Dim::toWolfram() const {
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);

	Vector<1> thisPos(this->getPosition());
	out << std::fixed;
	out << "GeometricTransformation[Disk[{0, 0}, {" << this->a << ", " << this->b << "}]," << std::endl;
	out << "    {RotationMatrix[" << this->getAngle() << "], {" << thisPos[0] << ", 0}}]";

	return out.str();
}

std::string Ellipse1Dim::toString() const {
    std::stringstream out;

	out.precision(std::numeric_limits< double >::max_digits10);
    Vector<1> thisPos(this->getPosition());
	out << "Ellipse{position: " << thisPos << "; a: " << this->a;
    out << "; b: " << this->b << "; angle: " << this->getAngle() << "}";
    return out.str();
}

std::string Ellipse1Dim::toPovray() const{
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);
	out << "  disc { <0.0, 0.0, 0.0002>, <0.0, 0.0, 1.0>, 1.0" << std::endl;;
	out << "	scale <" << this->a << ", " << this->b << ", 1.0>" << std::endl;
	out << "	rotate <0, 0, " << (180*this->getAngle()/M_PI) << ">" << std::endl;
	out << "	translate <";

	Vector<1> position = this->getPosition();
	for(unsigned short i=0; i<1; i++)
		out << position[i] << ", ";
	out << "0.0>" << std::endl << "	texture { pigment { color Red } }" << std::endl << "}" << std::endl;

	return out.str();
}

void Ellipse1Dim::store(std::ostream &f) const{
	Shape::store(f);
	f.write((char *)(&this->a), sizeof(double));
	f.write((char *)(&this->b), sizeof(double));
}

void Ellipse1Dim::restore(std::istream &f){
	Shape::restore(f);
	f.read((char *)(&this->a), sizeof(double));
	f.read((char *)(&this->b), sizeof(double));
	this->calculateU();

}

Shape<1, 1> *Ellipse1Dim::clone() const {
    return new Ellipse1Dim(*this);
}

bool Ellipse1Dim::pointInside(BoundaryConditions<1> *bc, const Vector<1> &position,
									 const Orientation<1> &orientation, const Orientation<1> &orientationRange) const {
	return this->pointInside(bc, position, orientation[0], orientation[0]+orientationRange[0]);
}

Matrix<2, 2> Ellipse1Dim::getRotationMatrix() const {
	return Matrix<2, 2>::rotation(this->getAngle());
}

Matrix<2, 2> Ellipse1Dim::getAntiRotationMatrix() const {
	return Matrix<2, 2>::rotation(-this->getAngle());
}

void Ellipse1Dim::normalizeAngleRange(double *angleFrom, double *angleTo, double interval) const {
	while (*angleFrom < this->getAngle()) {
		*angleFrom += interval;
		*angleTo += interval;
	}
	while (*angleFrom > this->getAngle() + interval) {
		*angleFrom -= interval;
		*angleTo -= interval;
	}
}

double Ellipse1Dim::normalizeAngle(double angle, double interval) const {
	while (angle < 0)
		angle += interval;
	while (angle > interval)
		angle -= interval;
	return angle;
}

void Ellipse1Dim::setOrientation(const Orientation<1> &orientation) {
	this->setAngle(orientation[0]);
}

double Ellipse1Dim::getAngle() const{
	return this->getOrientation()[0];
}
