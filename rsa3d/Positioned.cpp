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
}

Positioned::~Positioned() {
	delete this->position;
}

double* Positioned::getPosition(){
	return this->position;

}
