/*
 * Positioned.h
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#ifndef POSITIONED_H_
#define POSITIONED_H_

class Positioned {
protected:
	double *position;
	unsigned char dimension;

public:
	Positioned(unsigned char dimension);
	Positioned(const Positioned & other);
	virtual ~Positioned();

	virtual Positioned & operator=(const Positioned & other);

	// returns position of the shape's center
	double* getPosition();

	double distanceOf(Positioned & pos);
};

#endif /* POSITIONED_H_ */
