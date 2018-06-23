/*
 * Ellipse.cpp
 *
 *  Created on: 10.08.2017
 *      Author: ciesla
 */

#include "Ellipse.h"
#include "../../Utils.h"
#include "../../Intersection.h"


static const double EPSILON = 0.0000000001;

double Ellipse::longSemiAxis;
double Ellipse::shortSemiAxis;

void Ellipse::calculateU() {
	this->u[0]  = cos(this->getAngle()); this->u[1]  = sin(this->getAngle());
	this->uT[0] = -sin(this->getAngle()); this->uT[1] =  cos(this->getAngle());
}

Ellipse::Ellipse() : AnisotropicShape2D(){
	this->a = Ellipse::longSemiAxis;
	this->b = Ellipse::shortSemiAxis;
	this->setAngle(0);
}

void Ellipse::initClass(const std::string &args){
	double ratio = std::stod(args);
	Ellipse::shortSemiAxis = sqrt(1.0/(M_PI*ratio));
	Ellipse::longSemiAxis = ratio*shortSemiAxis;
	Shape<2,1>::setNeighbourListCellSize(2*longSemiAxis);
	Shape<2,1>::setVoxelSpatialSize(1.4*shortSemiAxis);
	Shape<2,1>::setVoxelAngularSize(M_PI);
	Shape<2,1>::setDefaultCreateShapeImpl <Ellipse> ();
}

double Ellipse::calculateF(double* r, double g) const {
	double d1 = (r[0]*this->u[0] + r[1]*this->u[1])  /this->a;
	double d2 = (r[0]*this->uT[0] + r[1]*this->uT[1]) /this->b;

	return 1 + g - d1*d1 - d2*d2;
}

void Ellipse::setAngle(double angle){
	AnisotropicShape2D::setAngle(angle);
	this->calculateU();
}

bool Ellipse::overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const {
	Ellipse es = dynamic_cast<const Ellipse&>(*s);
    this->applyBC(bc, &es);
	double d;
	d = (this->a/this->b - this->b/this->a)*sin(this->getAngle() - es.getAngle());
	double g = 2 + d*d;
	double r[] = {this->getPosition()[0] - es.getPosition()[0], this->getPosition()[1] - es.getPosition()[1]};
	double f1 = this->calculateF(r, g);
	double f2 = es.calculateF(r, g);
	double psi = 4*(f1*f1-3*f2)*(f2*f2-3*f1) - (9-f1*f2)*(9-f1*f2);
	if (psi>0 && (f1<0 || f2<0)){
		return false;
	}
	return true;
}

double Ellipse::getVolume() const {
	return M_PI*this->a*this->b;
}

// TODO not working properly

/*bool Ellipse::pointInside(BoundaryConditions *bc, double* da) const {
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

bool Ellipse::pointInside(BoundaryConditions<2> *bc, const Vector<2> &other, double angleFrom, double angleTo) const
{
    // Check exlusion zones for angle interval endpoinds - angleFrom and angleTo
    Ellipse ellTmp;
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

bool Ellipse::pointInsideSpecialArea(BoundaryConditions<2> *bc, const Vector<2> &other, double angleFrom,
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
    Vector<2> translationVector = bc->getTranslation(this->getPosition(), other);
    Vector<2> thisPos(this->getPosition());
    Vector<2> otherPos = Vector<2>(other) + translationVector;
	Vector<2> otherAligned = getAntiRotationMatrix() * (thisPos - otherPos);
	if (angleTo < this->getAngle())
		return pointInsideUnrotated(otherAligned, angleFrom - this->getAngle() + M_PI, angleTo - this->getAngle() + M_PI);
	else if (angleFrom < this->getAngle())
		return pointInsideUnrotated(otherAligned, angleFrom - this->getAngle() + M_PI, M_PI) &&
               pointInsideUnrotated(otherAligned, 0, angleTo - this->getAngle());
	else    // angleTo >= angle && angleFrom >= angle
		return pointInsideUnrotated(otherAligned, angleFrom - this->getAngle(), angleTo - this->getAngle());
}

bool Ellipse::pointInsideUnrotated(const Vector<2> &p, double angleFrom, double angleTo) const
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

bool Ellipse::withinAngle(const Vector<2> &p, double angleFrom, double angleTo) const
{
	double minT = atan2(b * sin(angleFrom), a * cos(angleFrom));
	double maxT = atan2(b * sin(angleTo), a * cos(angleTo));

	Vector<2> edgeMin{{cos(minT) * a, sin(minT) * b}};
	Vector<2> edgeMax{{cos(maxT) * a, sin(maxT) * b}};
	Vector<2> edgeCross = intersection::line_line(edgeMin, angleFrom, edgeMax, angleTo);
	double zoneAngle = getAngleToOrigin(p - edgeCross);
    
    return zoneAngle >= angleFrom && zoneAngle <= angleTo;
}

bool Ellipse::withinAngleCheckCollision(const Vector<2> &p, double lowerAngle, double upperAngle) const
{
	double minT = atan2(b * sin(lowerAngle), a * cos(lowerAngle));
	double maxT = atan2(b * sin(upperAngle), a * cos(upperAngle));

	return circleCollision(p, minT, maxT);
}

bool Ellipse::circleCollision(const Vector<2> &p, double tMin, double tMax) const
{
    Vector<2> start{{cos(tMin) * a, sin(tMin) * b}};
    Vector<2> end{{cos(tMax) * a, sin(tMax) * b}};
    if (intersection::line_circle(p, b, start, end))
        return true;

    double slope1 = atan2(a * sin(tMin), b * cos(tMin)) + M_PI/2;
    double slope2 = atan2(a * sin(tMax), b * cos(tMax)) + M_PI/2;
    Vector<2> intersect = intersection::line_line(start, slope1, end, slope2);
	if (!intersection::line_circle(p, b, start, intersect) && !intersection::line_circle(p, b, intersect, end))
		return false;

	return circleCollision(p, tMin, (tMin + tMax) / 2) || circleCollision(p, (tMin + tMax) / 2, tMax);
}

std::string Ellipse::toWolfram() const {
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);

	Vector<2> thisPos(this->getPosition());
	out << std::fixed;
	out << "GeometricTransformation[Disk[{0, 0}, {" << this->a << ", " << this->b << "}]," << std::endl;
	out << "    {RotationMatrix[" << this->getAngle() << "], " << thisPos << "}]";

	return out.str();
}

std::string Ellipse::toString() const {
    std::stringstream out;

	out.precision(std::numeric_limits< double >::max_digits10);
    Vector<2> thisPos(this->getPosition());
	out << "Ellipse{position: " << thisPos << "; a: " << this->a;
    out << "; b: " << this->b << "; angle: " << this->getAngle() << "}";
    return out.str();
}

std::string Ellipse::toPovray() const{
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);
	out << "  disc { <0.0, 0.0, 0.0002>, <0.0, 0.0, 1.0>, 1.0" << std::endl;;
	out << "	scale <" << this->a << ", " << this->b << ", 1.0>" << std::endl;
	out << "	rotate <0, 0, " << (180*this->getAngle()/M_PI) << ">" << std::endl;
	out << "	translate <";

	Vector<2> position = this->getPosition();
	for(unsigned short i=0; i<2; i++)
		out << position[i] << ", ";
	out << "0.0>" << std::endl << "	texture { pigment { color Red } }" << std::endl << "}" << std::endl;

	return out.str();
}

void Ellipse::store(std::ostream &f) const{
	Shape::store(f);
	f.write((char *)(&this->a), sizeof(double));
	f.write((char *)(&this->b), sizeof(double));
}

void Ellipse::restore(std::istream &f){
	Shape::restore(f);
	f.read((char *)(&this->a), sizeof(double));
	f.read((char *)(&this->b), sizeof(double));
	this->calculateU();

}

Shape<2, 1> *Ellipse::clone() const {
    return new Ellipse(*this);
}
