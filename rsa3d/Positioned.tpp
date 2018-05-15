/*
 * Positioned.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#include<algorithm>


template <unsigned short SPATIAL_DIMENSION>
int Positioned<SPATIAL_DIMENSION>::offset[(1 << SPATIAL_DIMENSION)][SPATIAL_DIMENSION];

template <unsigned short SPATIAL_DIMENSION>
double* Positioned<SPATIAL_DIMENSION>::getPosition() const {
    // TODO const_cast
	return (const_cast<Positioned*>(this))->position.data();
}

template<unsigned short SPATIAL_DIMENSION>
void Positioned<SPATIAL_DIMENSION>::setPosition(const double *position) {
    std::copy(position, position + SPATIAL_DIMENSION, this->position.begin());
}

template<unsigned short SPATIAL_DIMENSION>
void Positioned<SPATIAL_DIMENSION>::translate(double *v){
	double position[SPATIAL_DIMENSION];
	for(unsigned short i=0; i<SPATIAL_DIMENSION; i++)
		position[i] = this->getPosition()[i] + v[i];
	this->setPosition(position);
}
