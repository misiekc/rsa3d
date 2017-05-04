/*
 * Voxel.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#include <string.h>
#include "Voxel.h"
#include "Positioned.h"
#include <algorithm>

Voxel::Voxel(unsigned char dim) : Positioned(dim){
	this->index = 0;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
}

Voxel::Voxel(unsigned char dim, double* da, double s, int i) : Positioned(dim){
	std::copy(da, da+this->dimension, this->position);
	this->index = i;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
}

Voxel::Voxel(const Voxel & other) : Positioned(other){
	this->index = other.index;
	this->missCounter = other.missCounter;
	this->lastAnalyzed = other.lastAnalyzed;
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

