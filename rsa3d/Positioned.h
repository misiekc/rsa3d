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
private:
    double position[SPATIAL_DIMENSION];

protected:
    virtual void setPosition(const double *position);

public:
	Positioned();
	virtual ~Positioned();

	virtual Positioned<SPATIAL_DIMENSION> & operator=(const Positioned<SPATIAL_DIMENSION> & other);

	// returns position of the shape's center
	double* getPosition() const;

	//	Vector<SPATIAL_DIMENSION> getVectorPosition() const;
};

#include "Positioned.tpp"

#endif /* POSITIONED_H_ */
