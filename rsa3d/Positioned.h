/*
 * Positioned.h
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#ifndef POSITIONED_H_
#define POSITIONED_H_

// #include "Vector.h"

template <unsigned short SPATIAL_DIMENSION>
class Positioned {

protected:
	double position[SPATIAL_DIMENSION];

public:
	Positioned();
	Positioned(const Positioned<SPATIAL_DIMENSION> & other);
	virtual ~Positioned();

	virtual Positioned<SPATIAL_DIMENSION> & operator=(const Positioned<SPATIAL_DIMENSION> & other);

	// returns position of the shape's center
	double* getPosition();

	//	Vector<SPATIAL_DIMENSION> getVectorPosition() const;
};

#include "Positioned.tpp"

#endif /* POSITIONED_H_ */
