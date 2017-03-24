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
	int dimension;

public:
	Positioned(int dimension);
	virtual ~Positioned();

	// returns position of the shape's center
	double* getPosition();
};

#endif /* POSITIONED_H_ */
