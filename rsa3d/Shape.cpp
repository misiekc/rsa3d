/*
 * Shape.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include "Shape.h"
#include "Positioned.h"
#include <iostream>

Shape::Shape(const unsigned short dimension) : Positioned(dimension){
	this->no = 0;
	this->time = 0.0;
}

Shape::Shape(const Shape & other) : Positioned(other){
	this->no = 0;
	this->time = 0.0;
}

Shape::~Shape() {

}

void Shape::translate(double* v){
	for(unsigned short i=0; i<this->dimension; i++){
		this->position[i] += v[i];
	}
}

double Shape::minDistance(Shape *s){
	return 0.0;
}

std::string Shape::toString(){
	return "";
}

std::string Shape::toPovray() const{
	return "";
}

void Shape::store(std::ostream &f) const{
	f.write((char *)(&this->dimension), sizeof(unsigned char));
	f.write((char *)(&this->no), sizeof(int));
	f.write((char *)(&this->time), sizeof(double));
	f.write((char *)(this->position), this->dimension*sizeof(double));
}

void Shape::restore(std::istream &f){
	unsigned char d;
	f.read((char *)(&d), sizeof(unsigned char));

	if (f.gcount()==0){ // end of file
		return;
	}
	if (this->dimension!=d){
		std::cout << "[ERROR] cannot restore: incompatible dimensions: read " << f.gcount() << " bytes." << std::endl;
		return;
	}
	f.read((char *)(&this->no), sizeof(int));
	f.read((char *)(&this->time), sizeof(double));
	f.read((char *)(this->position), d*sizeof(double));
}


