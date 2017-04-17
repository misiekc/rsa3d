/*
 * Sphere.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#include <math.h>
#include <string>

#include "Sphere.h"
#include "../BoundaryConditions.h"

double Sphere::radius;
double Sphere::neighbourListCellSize;
double Sphere::voxelSize;
int Sphere::dimension;


Sphere::~Sphere() {
}


double Sphere::gamma(int d){
	double result;
	if(d%2==0){
		result = Sphere::g20;
		for (int i=4; i<=d; i+=2){
			result = (i/2.0)*result;
		}
	}else{ // d%2==1
		result = Sphere::g15;
		for (int i=3; i<=d; i+=2){
			result = (i/2.0)*result;
		}
	}
	return result;
}

double Sphere::volume(int d){
	return pow(M_PI, d/2.0) / gamma(d);
}

void Sphere::initClass(const std::string &args){
	Sphere::dimension = std::stoi(args);
	Sphere::radius = pow(1.0/Sphere::volume(Sphere::dimension), 1.0/Sphere::dimension);
	Sphere::neighbourListCellSize = 2.0*Sphere::radius;
	Sphere::voxelSize = pow(1.96, 1.0/Sphere::dimension)*Sphere::radius;
}

Shape * Sphere::create(RND *rnd){
	return new Sphere();
}

Sphere::Sphere() : Shape(Sphere::dimension){
	this->r = Sphere::radius;
}


double Sphere::getNeighbourListCellSize() {
	return Sphere::neighbourListCellSize;
}

double Sphere::getVoxelSize() {
	return Sphere::voxelSize;
}

int Sphere::overlap(BoundaryConditions *bc, Shape *s) {
	Sphere *sd = (Sphere *) s;
	double d2 = bc->distance2(this->position, sd->position);
	double r2 = this->r + sd->r;
	r2 *= r2;

	return (d2 < r2);
}

double Sphere::getVolume() {
	return Sphere::volume(Sphere::dimension)*pow(this->r, Sphere::dimension);
}

int Sphere::pointInside(BoundaryConditions *bc, double* da) {
	double d2;
	if (bc!=NULL)
		d2 = bc->distance2(da, this->position);
	else{
		d2 = 0.0;
		for(int i=0; i<Sphere::dimension; i++)
			d2 += (da[i]-this->position[i])*(da[i]-this->position[i]);
	}
	return (d2<4.0*this->r*this->r);
}

