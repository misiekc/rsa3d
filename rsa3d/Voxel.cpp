/*
 * Voxel.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#include <string.h>
#include "Voxel.h"

Voxel::Voxel(int dim, double* da, double s, int i){
	this->dimension = dim;
	this->center = new double[this->dimension];
	memcpy(this->center, da, sizeof(double)*this->dimension);
	this->index = i;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
	}

Voxel::~Voxel() {
	delete this->dimension;
}

void Voxel::miss(){
	this->missCounter++;
}

int Voxel::getMissCounter(){
	return this->missCounter;
}

void Voxel::resetMissCounter(){
	this->missCounter = 0;
}

double* Voxel::getPosition(){
	return this->center;
}

