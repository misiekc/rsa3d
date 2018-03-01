/*
 * Positioned.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#include<algorithm>

template <unsigned short SPATIAL_DIMENSION>
Positioned<SPATIAL_DIMENSION>::Positioned() {
	for(unsigned short i=0; i<SPATIAL_DIMENSION; i++){
		this->position[i] = 0.0;
	}
}

template <unsigned short SPATIAL_DIMENSION>
Positioned<SPATIAL_DIMENSION>::Positioned(const Positioned & other) {
    std::copy(other.position, other.position+SPATIAL_DIMENSION, this->position);
}


template <unsigned short SPATIAL_DIMENSION>
Positioned<SPATIAL_DIMENSION>::~Positioned() {
}

template <unsigned short SPATIAL_DIMENSION>
Positioned<SPATIAL_DIMENSION> & Positioned<SPATIAL_DIMENSION>::operator=(const Positioned<SPATIAL_DIMENSION> & other){
    // Self assingment, skip
    if (this == &other)
        return *this;

    std::copy(other.position, other.position+SPATIAL_DIMENSION, this->position);
    return *this;
}

template <unsigned short SPATIAL_DIMENSION>
double* Positioned<SPATIAL_DIMENSION>::getPosition(){
	return this->position;
}

