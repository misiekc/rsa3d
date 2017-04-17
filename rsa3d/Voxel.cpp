/*
 * Voxel.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#include <string.h>
#include "Voxel.h"
#include "Positioned.h"

Voxel::Voxel(int dim) : Positioned(dim){
	this->index = 0;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
}

Voxel::Voxel(int dim, double* da, double s, int i) : Positioned(dim){
	this->dimension = dim;
	memcpy(this->position, da, sizeof(double)*this->dimension);
	this->index = i;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
}

Voxel::~Voxel() {
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

