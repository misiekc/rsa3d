/*
 * Polygon.cpp
 *
 *  Created on: 16.04.2018
 *      Author: ciesla
 */

#include "Polygon.h"
#include <cmath>
#include <climits>


double Polygon::findVoxelSpatialSize(Polygon *p){
	double *pos = p->getPosition();

	double d, dMin = p->segments[0]->distance(pos[0],  pos[1]);
	for(int i=1; i<p->n; i++){
		d = p->segments[i]->distance(pos[0],  pos[1]);
		if (d<dMin)
			dMin = d;
	}
	return dMin/sqrt(2.0);
}

double Polygon::findNeighbourGridSize(Polygon *p){
	double *pos = p->getPosition();

	double d, dx, dy, dMax = 0.0;
	dMax = p->segments[0]->distance(pos[0],  pos[1]);
	for(int i=1; i<p->n; i++){
		dx = pos[0] - p->segments[i]->getX1();
		dy = pos[1] - p->segments[i]->getY1();
		d = sqrt(dx*dx + dy*dy);
		if (d>dMax)
			dMax = d;
	}
	return 2*dMax;

}



Polygon::Polygon(unsigned short n) : Shape(){
	this->n = n;
	this->segments = new Segment2D*[n];
}

Polygon::Polygon(const Polygon &p){
	this->n = p.n;
	this->segments = new Segment2D*[n];
	std::copy(p.segments, p.segments + p.n, this->segments);
}


Polygon::~Polygon() {
	for(unsigned short i=0; i<this->n; i++){
		delete this->segments[i];
	}
	delete[] this->segments;
}

void Polygon::applyBC(BoundaryConditions *bc, Shape<2,0> *second) const {
	Polygon *p = (Polygon *)second;
    double translation[2];
    bc->getTranslation(translation, this->getPosition(), second->getPosition());
    second->translate(translation);
    for(int i=0; i<p->n; i++){
    	p->segments[i]->translate(translation);
    }
}


int Polygon::overlap(BoundaryConditions *bc, Shape<2, 0> *s) const{
	Polygon p( *((Polygon *)s));
	this->applyBC(bc, &p);

	for(int i=0; i<this->n; i++){
		for(int j=0; j<p.n; j++){
			if(this->segments[i]->isCrossing(p.segments[j]))
				return 1;
		}
	}
	return 0;
}

int Polygon::pointInside(BoundaryConditions *bc, double* position, const std::array<double, 0> &orientation, double orientationRange) const{
	double translation[2], pos[2];
	bc->getTranslation(translation,  this->getPosition(), position);
	pos[0] = position[0] + translation[0];
	pos[1] = position[1] + translation[1];

	double d, dMin = this->segments[0]->distance(pos[0],  pos[1]);
	for(int i=1; i<this->n; i++){
		d = this->segments[i]->distance(pos[0],  pos[1]);
		if (d<dMin)
			dMin = d;
	}

	if (dMin < 2*sqrt(2.0) * (Shape<2,0>::getVoxelSpatialSize()) )
		return 1;
	else
		return 0;
}
