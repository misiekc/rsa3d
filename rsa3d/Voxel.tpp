/*
 * Voxel.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

template <ushort DIMENSION>
Voxel<DIMENSION>::Voxel(){
	this->index = 0;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
}

template <ushort DIMENSION>
Voxel<DIMENSION>::Voxel(double* da, double s, int i){
	std::copy(da, da+DIMENSION, this->position);
	this->index = i;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
}

template <ushort DIMENSION>
Voxel<DIMENSION>::~Voxel() {
}


template <ushort DIMENSION>
Voxel<DIMENSION>::Voxel(const Voxel & other){
	std::copy(other.position, other.position+DIMENSION, this->position);
	this->index = other.index;
	this->missCounter = other.missCounter;
	this->lastAnalyzed = other.lastAnalyzed;
}

template <ushort DIMENSION>
void Voxel<DIMENSION>::miss(){
	this->missCounter++;
}

template <ushort DIMENSION>
int Voxel<DIMENSION>::getMissCounter(){
	return this->missCounter;
}

template <ushort DIMENSION>
void Voxel<DIMENSION>::resetMissCounter(){
	this->missCounter = 0;
}

template <ushort DIMENSION>
bool Voxel<DIMENSION>::isInside(double *da, double size){
	for(int i=0; i<DIMENSION; i++){
		if (da[i]<this->position[i])
			return false;
		if (da[i]>=(this->position[i]+size))
			return false;
	}
	return true;
}

template <ushort DIMENSION>
double* Voxel<DIMENSION>::getPosition(){
	return this->position;
}
