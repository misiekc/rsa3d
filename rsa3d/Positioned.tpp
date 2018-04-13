/*
 * Positioned.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#include<algorithm>


template <unsigned short SPATIAL_DIMENSION>
double* Positioned<SPATIAL_DIMENSION>::getPosition() const {
    // TODO const_cast
	return (const_cast<Positioned*>(this))->position.data();
}

template<unsigned short SPATIAL_DIMENSION>
void Positioned<SPATIAL_DIMENSION>::setPosition(const double *position) {
    std::copy(position, position + SPATIAL_DIMENSION, this->position.begin());
}

