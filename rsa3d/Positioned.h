/*
 * Positioned.h
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#ifndef POSITIONED_H_
#define POSITIONED_H_

template <unsigned short DIMENSION>
class Positioned {

protected:
	double position[DIMENSION];

public:
	Positioned();
	Positioned(const Positioned<DIMENSION> & other);
	virtual ~Positioned();

	virtual Positioned<DIMENSION> & operator=(const Positioned<DIMENSION> & other);

	// returns position of the shape's center
	double* getPosition();
};

#include "Positioned.tpp"

#endif /* POSITIONED_H_ */
