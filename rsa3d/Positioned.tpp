/*
 * Positioned.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#include<algorithm>

template <unsigned short DIMENSION>
Positioned<DIMENSION>::Positioned() {
	for(unsigned short i=0; i<DIMENSION; i++){
		this->position[i] = 0.0;
	}
}

template <unsigned short DIMENSION>
Positioned<DIMENSION>::Positioned(const Positioned & other) {
    std::copy(other.position, other.position+DIMENSION, this->position);
}


template <unsigned short DIMENSION>
Positioned<DIMENSION>::~Positioned() {
}

template <unsigned short DIMENSION>
Positioned<DIMENSION> & Positioned<DIMENSION>::operator=(const Positioned<DIMENSION> & other){
    // Self assingment, skip
    if (this == &other)
        return *this;

    std::copy(other.position, other.position+DIMENSION, this->position);
    return *this;
}

template <unsigned short DIMENSION>
double* Positioned<DIMENSION>::getPosition(){
	return this->position;
}
