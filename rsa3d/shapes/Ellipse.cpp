/*
 * Ellipse.cpp
 *
 *  Created on: 10.08.2017
 *      Author: ciesla
 */

#include "Ellipse.h"
#include "../Utils.h"
#include "../tests/utility/MockBC.h"


#define EPSILON 0.0000001

double Ellipse::longSemiAxis;
double Ellipse::shortSemiAxis;
double Ellipse::neighbourListCellSize;
double Ellipse::voxelSize;

void Ellipse::calculateU(){
	this->u[0]  = cos(this->angle); this->u[1]  = sin(this->angle);
	this->uT[0] = -sin(this->angle); this->uT[1] =  cos(this->angle);
}

Ellipse::Ellipse() : AnisotropicShape2D(){
	this->a = Ellipse::longSemiAxis;
	this->b = Ellipse::shortSemiAxis;
	this->angle = 0.0;
	this->calculateU();
}

Ellipse::Ellipse(const Ellipse &other) : AnisotropicShape2D(other), a(other.a), b(other.b) {
    for(unsigned char i = 0; i < 2; i++){
        this->u[i]  = other.u[i];
        this->uT[i] = other.uT[i];
    }
}

Ellipse & Ellipse::operator = (const Ellipse & el){
	this->a = el.a;
	this->b = el.b;
	this->angle = el.angle;
	for(unsigned char i=0; i<2; i++){
		this->u[i]  = el.u[i];
		this->uT[i] = el.uT[i];
	}
	return *this;
}

void Ellipse::initClass(const std::string &args){
	double ratio = std::stod(args);
	Ellipse::shortSemiAxis = sqrt(1.0/(M_PI*ratio));
	Ellipse::longSemiAxis = ratio*shortSemiAxis;
	Ellipse::neighbourListCellSize = 2*longSemiAxis;
	Ellipse::voxelSize = 1.4*shortSemiAxis;
}

Shape<2> * Ellipse::create(RND *rnd){
	Ellipse *el = new Ellipse();
	el->a = Ellipse::longSemiAxis;
	el->b = Ellipse::shortSemiAxis;
	el->angle = rnd->nextValue()*2*M_PI;
	el->calculateU();
	return el;
}

double Ellipse::calculateF(double* r, double g){
	double d1 = (r[0]*this->u[0] + r[1]*this->u[1])  /this->a;
	double d2 = (r[0]*this->uT[0] + r[1]*this->uT[1]) /this->b;

	return 1 + g - d1*d1 - d2*d2;
}

int Ellipse::overlap(BoundaryConditions *bc, Shape *s) {
	Ellipse es = *((Ellipse *)s);
	double da[2];
	bc->getTranslation(da, this->position, es.position);
	es.translate(da);
	double d;
	d = (this->a/this->b - this->b/this->a)*sin(this->angle - es.angle);
	double g = 2 + d*d;
	double r[] = {this->position[0] - es.position[0], this->position[1] - es.position[1]};
	double f1 = this->calculateF(r, g);
	double f2 = es.calculateF(r, g);
	double psi = 4*(f1*f1-3*f2)*(f2*f2-3*f1) - (9-f1*f2)*(9-f1*f2);
	if (psi>0 && (f1<0 || f2<0)){
		return false;
	}
	return true;
}

double Ellipse::getVolume() {
	return M_PI*this->a*this->b;
}


int Ellipse::pointInside(BoundaryConditions *bc, double* da) {
	double ta[2];
	double tmp[2];

	tmp[0] = da[0]; tmp[1] = da[1];
	bc->getTranslation(ta, this->position, tmp);
	da[0] += ta[0];
	da[1] += ta[1];

	da[0] -= this->position[0];
	da[1] -= this->position[1];

	rotate(da, -this->angle);

	double dx = da[0]/(this->a+this->b);
	double dy = da[1]/(2*this->b);
	return (dx*dx+dy*dy < 1);
}

bool lineCircleCollision(const Vector<2> &center, double r, const Vector<2> &p1, const Vector<2> &p2)
{
    double a = p2[0] - p1[0];
    double b = p2[1] - p1[1];
    double c = center[0] - p1[0];
    double d = center[1] - p1[1];
    if ( (d*a - c*b)*(d*a - c*b) <= r*r * (a*a + b*b) ) {
        if (c*c + d*d <= r*r) //first end
            return true;
        else if ( (a-c)*(a-c) + (b-d)*(b-d) <= r*r) //second end
            return true;
        else if ( (c*a + d*b >= 0) && (c*a + d*b <= a*a + b*b))
            return true;
        else
            return false;
    } else {
        return false;
    }
}

Vector<2> lineLineIntersection(const Vector<2> &line1, double line1Angle, const Vector<2> &line2, double line2Angle)
{
    double a = tan(line1Angle);
    double b = tan(line2Angle);
    double c = line1[1] - a * line1[0];
    double d = line2[1] - b * line2[0];

    double a_b = a - b;
    if( a_b == 0 )//parallel lines
        return Vector<2>{{line1[0], line1[1]}};
    else
        return Vector<2>{{(d - c) / a_b, (a * d - b * c) / a_b}};
}

int Ellipse::pointInside(BoundaryConditions *bc, double *other, double angleFrom, double angleTo)
{
    // to be inside an exclusion zone the point have to bee inside all exclusion zones between angleFrom and angleTo
    Ellipse ellTmp;
    ellTmp.translate(other);

    // checking if the point is inside the exclusion zone for angleFrom
    ellTmp.setAngle(angleFrom);
    if (!this->overlap(bc, &ellTmp))
        return false;

    // checking if the point is inside the exclusion zone for angleTo
    ellTmp.setAngle(angleTo);
    if (!this->overlap(bc, &ellTmp))
        return false;

	return pointInside0(bc, other, angleFrom, angleTo);
}

int Ellipse::pointInside0(BoundaryConditions *bc, double *other, double angleFrom, double angleTo) {
	this->normalizeAngleRange(angleFrom, angleTo, M_PI);

    // If angleTo not in normal range (see normalizeAngleRange), divide it and check separately
    if (angleTo > this->angle + M_PI) {
        return pointInside0(bc, other, angleFrom, this->angle + M_PI - EPSILON) &&
               pointInside0(bc, other, this->angle, angleTo - M_PI);
    }
	
	Vector<2> p = getAntiRotationMatrix() * (this->getVectorPosition() - this->applyBC(bc, other));
	if (angleTo < angle)
		return withinExclusionZoneUnrotated(p, angleFrom - angle + M_PI, angleTo - angle + M_PI);
	else if (angleFrom < angle)
		return withinExclusionZoneUnrotated(p, angleFrom - angle + M_PI, M_PI) &&
			   withinExclusionZoneUnrotated(p, 0, angleTo - angle);
	else
		return withinExclusionZoneUnrotated(p, angleFrom - angle, angleTo - angle);
}

bool Ellipse::withinExclusionZoneUnrotated(const Vector<2> &p, double lowerAngle, double upperAngle) const
{
	if (lowerAngle <= M_PI/2 && upperAngle >= M_PI/2) 
    {
		if (p[0] <= 0 && p[1] >= 0 && withinAngle(p, lowerAngle + M_PI/2, M_PI))
			return withinAngleCheckCollision(p, lowerAngle + M_PI/2, M_PI);
		if (p[0] <= 0 && p[1] <= 0 && withinAngle(p, M_PI + EPSILON, upperAngle + M_PI/2))
			return withinAngleCheckCollision(p, M_PI + EPSILON, upperAngle + M_PI/2);
		if (p[0] >= 0 && p[1] <= 0 && withinAngle(-p, lowerAngle + M_PI/2, M_PI))
			return withinAngleCheckCollision(-p, lowerAngle + M_PI/2, M_PI);
		if (p[0] >= 0 && p[1] >= 0 && withinAngle(-p, M_PI + EPSILON, upperAngle + M_PI/2))
			return withinAngleCheckCollision(-p, M_PI + EPSILON, upperAngle + M_PI/2);
	} 
    else if (lowerAngle <= M_PI/2 && upperAngle <= M_PI/2) 
    {
		if (p[0] <= 0 && p[1] >= 0 && withinAngle(p, lowerAngle + M_PI/2, upperAngle + M_PI/2))
			return withinAngleCheckCollision(p, lowerAngle + M_PI/2, upperAngle + M_PI/2);
		if (p[0] >= 0 && p[1] <= 0 && withinAngle(-p, lowerAngle + M_PI/2, upperAngle + M_PI/2))
			return withinAngleCheckCollision(-p, lowerAngle + M_PI/2, upperAngle + M_PI/2);
	} 
    else 
    {
		if (p[0] <= 0 && p[1] <= 0 && withinAngle(p, lowerAngle + M_PI/2, upperAngle + M_PI/2))
			return withinAngleCheckCollision(p, lowerAngle + M_PI/2, upperAngle + M_PI/2);
		if (p[0] >= 0 && p[1] >= 0 && withinAngle(-p, lowerAngle + M_PI/2, upperAngle + M_PI/2))
			return withinAngleCheckCollision(-p, lowerAngle + M_PI/2, upperAngle + M_PI/2);
	}

	return true;
}

bool Ellipse::withinAngle(const Vector<2> &p, double lowerAngle, double upperAngle) const
{
	double minT = atan2(b * sin(lowerAngle), a * cos(lowerAngle));
	double maxT = atan2(b * sin(upperAngle), a * cos(upperAngle));

	Vector<2> edgeMin{{cos(minT) * a, sin(minT) * b}};
	Vector<2> edgeMax{{cos(maxT) * a, sin(maxT) * b}};
	Vector<2> edgeCross = lineLineIntersection(edgeMin, lowerAngle, edgeMax, upperAngle);
	double zoneAngle = getAngleToOrigin(p - edgeCross);
    
    return zoneAngle >= lowerAngle && zoneAngle <= upperAngle;
}

bool Ellipse::withinAngleCheckCollision(const Vector<2> &p, double lowerAngle, double upperAngle) const
{
	double minT = atan2(b * sin(lowerAngle), a * cos(lowerAngle));
	double maxT = atan2(b * sin(upperAngle), a * cos(upperAngle));

	return testCircleEllipseCollision(p, minT, maxT);
}

bool Ellipse::testCircleEllipseCollision(const Vector<2> &p, double tMin, double tMax) const
{
    Vector<2> start{{cos(tMin) * a, sin(tMin) * b}};
    Vector<2> end{{cos(tMax) * a, sin(tMax) * b}};
    if (lineCircleCollision(p, b, start, end))
        return true;

    double kont1 = atan2(a * sin(tMin), b * cos(tMin)) + M_PI / 2;
    double kont2 = atan2(a * sin(tMax), b * cos(tMax)) + M_PI / 2;
    Vector<2>intersect = lineLineIntersection(start, kont1, end, kont2);
	if (!lineCircleCollision(p, b, start, intersect) && !lineCircleCollision(p, b, intersect, end))
		return false;
    
	return testCircleEllipseCollision(p, tMin, (tMin + tMax) / 2) || testCircleEllipseCollision(p, (tMin + tMax) / 2, tMax);
}

void Ellipse::setAngle(double d){
	this->angle = d;
	this->calculateU();
}

std::string Ellipse::toWolfram() const {
	std::stringstream out;

	out << std::fixed;
	out << "GeometricTransformation[Disk[{0, 0}, {" << this->a << ", " << this->b << "}]," << std::endl;
	out << "    {RotationMatrix[" << this->angle << "], " << this->getVectorPosition() << "}]";

	return out.str();
}

double Ellipse::getNeighbourListCellSize() {
	return Ellipse::neighbourListCellSize;
}

double Ellipse::getVoxelSize() {
	return Ellipse::voxelSize;
}

std::string Ellipse::toString() {
    std::stringstream out;
    out << "Ellipse{position: " << this->getVectorPosition() << "; a: " << this->a;
    out << "; b: " << this->b << "; angle: " << this->angle << "}";
    return out.str();
}

std::string Ellipse::toPovray() const{
	//TODO
	return "";
}

void Ellipse::store(std::ostream &f) const{
	Shape::store(f);
	f.write((char *)(&this->a), sizeof(double));
	f.write((char *)(&this->b), sizeof(double));
	f.write((char *)(&this->angle), sizeof(double));
}

void Ellipse::restore(std::istream &f){
	Shape::restore(f);
	f.read((char *)(&this->a), sizeof(double));
	f.read((char *)(&this->b), sizeof(double));
	f.read((char *)(&this->angle), sizeof(double));
	this->calculateU();

}
