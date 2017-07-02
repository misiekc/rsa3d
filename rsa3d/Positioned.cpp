/*
 * Positioned.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#include "Positioned.h"
#include <algorithm>
#include <cmath>

Positioned::Positioned(unsigned char dim) {
	this->dimension = dim;
	this->position = new double[dim];
	for(unsigned char i=0; i<dim; i++){
		this->position[i] = 0.0;
	}
}

Positioned::Positioned(const Positioned & other) {
	this->dimension = other.dimension;
	this->position = new double[this->dimension];
    std::copy(other.position, other.position+other.dimension, this->position);
}


Positioned::~Positioned() {
	delete[] this->position;
}

Positioned & Positioned::operator=(const Positioned & other){
    // Self assingment, skip
    if (this == &other)
        return *this;

    this->dimension = other.dimension;
    delete [] this->position;
    std::copy(other.position, other.position+other.dimension, this->position);
    return *this;
}

double* Positioned::getPosition(){
	return this->position;
}

double Positioned::distanceOf(Positioned &pos){
	double d, res = 0.0;
	if (this->dimension != pos.dimension){
		return -1;
	}
	for(unsigned char i=0; i<this->dimension; i++){
		d = (this->position[i] - pos.position[i]);
		res += d*d;
	}
	return sqrt(res);
}
