/*
 * Voxel.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::Voxel(){
	this->index = 0;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::Voxel(double* pos, double *angle, int i){
	std::copy(pos, pos+SPATIAL_DIMENSION, this->position);
	std::copy(angle, angle+ANGULAR_DIMENSION, this->orientation);
	this->index = i;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::~Voxel() {
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::Voxel(const Voxel & other){
	std::copy(other.position, other.position+SPATIAL_DIMENSION, this->position);
	std::copy(other.orientation, other.orientation+ANGULAR_DIMENSION, this->position);
	this->index = other.index;
	this->missCounter = other.missCounter;
	this->lastAnalyzed = other.lastAnalyzed;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::miss(){
	this->missCounter++;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getMissCounter(){
	return this->missCounter;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::resetMissCounter(){
	this->missCounter = 0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::isInside(double *pos, double size){
	for(int i=0; i<SPATIAL_DIMENSION; i++){
		if (pos[i]<this->position[i])
			return false;
		if (pos[i]>=(this->position[i]+size))
			return false;
	}
	return true;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::isInside(double *pos, double size, double *angle, double asize){
	for(int i=0; i<SPATIAL_DIMENSION; i++){
		if (pos[i]<this->position[i])
			return false;
		if (pos[i]>=(this->position[i]+size))
			return false;
	}

	for(int i=0; i<ANGULAR_DIMENSION; i++){
		if (angle[i]<this->orientation[i])
			return false;
		if (angle[i]>=(this->orientation[i]+size))
			return false;
	}

	return true;
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double* Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getPosition(){
	return this->position;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double* Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getOrientation(){
	return this->orientation;
}
