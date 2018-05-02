/*
 * Segment.cpp
 *
 *  Created on: 29-10-2013
 *      Author: ciesla
 */

#include <cmath>
#include "Segment2D.h"
#include "../../BoundaryConditions.h"

double Segment2D::det(double x1, double y1, double x2, double y2, double x3, double y3){
	return x1*y2 + y1*x3 + x2*y3 - y1*x2 - y2*x3 - x1*y3;
}

int Segment2D::onSegment(double x1, double y1, double x2, double y2, double x3, double y3){
	if (   ((std::min(x1, x2) <= x3) && (x3 <= std::max(x1, x2)))
       	&& ((std::min(y1, y2) <= y3) && (y3 <= std::max(y1, y2))) ) return 1;
	return 0;

}

int Segment2D::isCrossing(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
	double d1 = Segment2D::det(x1, y1, x2, y2, x3, y3);
	double d2 = Segment2D::det(x1, y1, x2, y2, x4, y4);
	double d3 = Segment2D::det(x3, y3, x4, y4, x1, y1);
	double d4 = Segment2D::det(x3, y3, x4, y4, x2, y2);
	if (d1*d2<0 && d3*d4<0) return 1;
	if (d1==0 && Segment2D::onSegment(x1, y1, x2, y2, x3, y3)) return 1;
	if (d2==0 && Segment2D::onSegment(x1, y1, x2, y2, x4, y4)) return 1;
	if (d3==0 && Segment2D::onSegment(x3, y3, x4, y4, x1, y1)) return 1;
	if (d4==0 && Segment2D::onSegment(x3, y3, x4, y4, x2, y2)) return 1;
	return 0;
}

void Segment2D::calculateXY(){
	double *position = this->getPosition();
	this->x1 = position[0] - this->radius*cos(this->orientation);
	this->x2 = position[0] + this->radius*cos(this->orientation);
	this->y1 = position[0] - this->radius*sin(this->orientation);
	this->y2 = position[0] + this->radius*sin(this->orientation);
}

Segment2D::Segment2D(){
	this->radius = 0.0;
	this->orientation = 0.0;
	this->x1 = 0.0; this->x2 = 0.0; this->y1 = 0.0; this->y2 = 0.0;
}

Segment2D::Segment2D(double x1, double y1, double x2, double y2){
	double position[2] = {0.5*(x1+x2), 0.5*(y1+y2)};
	this->setPosition(position);
	double dx = (x2-x1), dy = (y2-y1);
	this->radius = 0.5*sqrt(dx*dx + dy*dy);
	this->orientation = atan2(dy, dx);
	this->x1 = x1;
	this->x2 = x2;
	this->y1 = y1;
	this->y2 = y2;
}

int Segment2D::isCrossing(Segment2D *s){
	return Segment2D::isCrossing(this->x1, this->y1, this->x2, this->y2, s->x1, s->y1, s->x2, s->y2);
}

void Segment2D::setPosition(const double *v){
	const double *oldPos = this->getPosition();
	Positioned::setPosition(v);
	const double *newPos = this->getPosition();
    double diff[] = {newPos[0] - oldPos[0], newPos[1] - oldPos[1]};

	this->x1 += diff[0]; this->y1 += diff[1];
	this->x2 += diff[0]; this->y2 += diff[1];
}

void Segment2D::rotate(double phi){
	double x1, x2, y1, y2;

	x1 = this->x1*cos(phi) - this->y1*sin(phi);
	y1 = this->x1*sin(phi) + this->y1*cos(phi);

	x2 = this->x2*cos(phi) - this->y2*sin(phi);
	y2 = this->x2*sin(phi) + this->y2*cos(phi);

	this->x1 = x1;
	this->y1 = y1;
	this->x2 = x2;
	this->y2 = y2;
}

double Segment2D::length(){
	double dx = this->x2 - this->x1;
	double dy = this->y2 - this->y1;
	return sqrt(dx*dx + dy*dy);
}

double Segment2D::distance(double x, double y){
	double dx1p = (x - this->x1);
	double dy1p = (y - this->y1);
	double d1p = sqrt(dx1p*dx1p + dy1p*dy1p);

	double dx2p = (x - this->x2);
	double dy2p = (y - this->y2);
	double d2p = sqrt(dx2p*dx2p + dy2p*dy2p);

	double dx = (this->x2 - this->x1);
	double dy = (this->y2 - this->y1);

	double cos1 = (dx1p*dx + dy1p*dy)/(sqrt(dx*dx + dy*dy)*d1p);
	double cos2 = (dx2p*dx + dy2p*dy)/(sqrt(dx*dx + dy*dy)*d2p);

	if (cos1*cos2<0) // (x, y) "between" ends of the segment
		return d1p*sqrt(1-cos1*cos1);
	else
		return std::min(d1p, d1p);
}

void Segment2D::applyBC(BoundaryConditions *bc, Segment2D *second) const {

    double translation[2];
    bc->getTranslation(translation, this->getPosition(), second->getPosition());
    second->translate(translation);
}


double Segment2D::getX1(){
	return this->x1;
}

double Segment2D::getY1(){
	return this->y1;
}

double Segment2D::getX2(){
	return this->x2;
}

double Segment2D::getY2(){
	return this->y2;
}
