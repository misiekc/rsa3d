/*
 * Voxel.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#include <string.h>
#include "Voxel.h"
#include "VoxelList.h"
#include "Positioned.h"
#include <algorithm>

Voxel::Voxel(){
	this->position = new double[VoxelList::dimension]();
	this->index = 0;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
}

Voxel::Voxel(double* da, double s, int i){
	this->position = new double[VoxelList::dimension];
	std::copy(da, da+VoxelList::dimension, this->position);
	this->index = i;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
}

Voxel::Voxel(const Voxel & other){
	this->position = new double[VoxelList::dimension];
	std::copy(other.position, other.position+VoxelList::dimension, this->position);
	this->index = other.index;
	this->missCounter = other.missCounter;
	this->lastAnalyzed = other.lastAnalyzed;
}

Voxel::~Voxel() {
	delete[] this->position;
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

bool Voxel::isInside(double *da, double size){
	for(int i=0; i<VoxelList::dimension; i++){
		if (da[i]<this->position[i])
			return false;
		if (da[i]>=(this->position[i]+size))
			return false;
	}
	return true;
}

double* Voxel::getPosition(){
	return this->position;
}
