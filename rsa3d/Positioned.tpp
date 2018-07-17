/*
 * Positioned.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#include<algorithm>
#include "Utils.h"


template <unsigned short SPATIAL_DIMENSION>
int Positioned<SPATIAL_DIMENSION>::offset[(1 << SPATIAL_DIMENSION)][SPATIAL_DIMENSION];

template <unsigned short SPATIAL_DIMENSION>
const Vector<SPATIAL_DIMENSION> &Positioned<SPATIAL_DIMENSION>::getPosition() const {
	return this->position;
}

template<unsigned short SPATIAL_DIMENSION>
void Positioned<SPATIAL_DIMENSION>::prepareOffset() {
	std::array<int, SPATIAL_DIMENSION> in{};
	in.fill(0);
	int index = 0;
	do {
		std::copy(in.begin(), in.end(), Positioned<RSA_SPATIAL_DIMENSION>::offset[index]);
		index++;
	} while(increment(in.data(), RSA_SPATIAL_DIMENSION, 1));
}

template<unsigned short SPATIAL_DIMENSION>
void Positioned<SPATIAL_DIMENSION>::setPosition(const Vector<SPATIAL_DIMENSION> &position) {
	this->position = position;
}

template<unsigned short SPATIAL_DIMENSION>
void Positioned<SPATIAL_DIMENSION>::translate(const Vector<SPATIAL_DIMENSION> &v){
	this->setPosition(this->position + v);
}
