/*
 * Sphere.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#include <math.h>

#include "Sphere.h"
#include "../BoundaryConditions.h"

#include <iostream>

double Sphere::radius;
double Sphere::neighbourListCellSize;
double Sphere::voxelSize;
unsigned char Sphere::staticDimension;


Sphere::~Sphere() {
}


double Sphere::gamma(unsigned char d){
	double result;
	if(d%2==0){
		result = Sphere::g20;
		for (int i=4; i<=d; i+=2){
			result = (i/2.0)*result;
		}
	}else{ // d%2==1
		result = Sphere::g15;
		for (unsigned char i=3; i<=d; i+=2){
			result = (i/2.0)*result;
		}
	}
	return result;
}

double Sphere::volume(unsigned char d){
	return pow(M_PI, d/2.0) / gamma(d);
}

void Sphere::initClass(const std::string &args){
	Sphere::staticDimension = std::stoi(args);
	Sphere::radius = pow(1.0/Sphere::volume(Sphere::staticDimension), 1.0/Sphere::staticDimension);
	Sphere::neighbourListCellSize = 2.0*Sphere::radius;
	Sphere::voxelSize = pow(1.96, 1.0/Sphere::staticDimension)*Sphere::radius;
}

Shape * Sphere::create(RND *rnd){
	return new Sphere();
}

Sphere::Sphere() : Shape(Sphere::staticDimension){
	this->r = Sphere::radius;
}

Sphere::Sphere(unsigned char dim) : Shape(dim){
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
	return Sphere::volume(Sphere::staticDimension)*pow(this->r, Sphere::staticDimension);
}

int Sphere::pointInside(BoundaryConditions *bc, double* da) {
	double d2;
	if (bc!=NULL)
		d2 = bc->distance2(da, this->position);
	else{
		d2 = 0.0;
		for(unsigned char i=0; i<this->dimension; i++)
			d2 += (da[i]-this->position[i])*(da[i]-this->position[i]);
	}
	return (d2<4.0*this->r*this->r);
}

std::string Sphere::toPovray(){
	std::string s = "  sphere { < ";
	for(unsigned char i=0; i<this->dimension; i++)
		s += std::to_string(this->position[i]) + ", ";
	s += "0.0>, " + std::to_string(this->r) +"\n    texture { pigment { color Red } }\n  }\n";
	return s;
}

void Sphere::store(std::ostream &f){
	Shape::store(f);
	f.write((char *)(&this->r), sizeof(double));
}

void Sphere::restore(std::istream &f){
	Shape::restore(f);
	f.read((char *)(&this->r), sizeof(double));
}

