/*
 * Positioned.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#include "Positioned.h"

Positioned::Positioned(int dim) {
	this->dimension = dim;
	this->position = new double[dim];
	for(int i=0; i<dim; i++){
		this->position[i] = 0.0;
	}
}

Positioned::~Positioned() {
	delete[] this->position;
}

double* Positioned::getPosition(){
	return this->position;

}
