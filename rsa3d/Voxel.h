/*
 * Voxel.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef VOXEL_H_
#define VOXEL_H_

#include "Positioned.h"

class Voxel{

private:
	double *position;
	short missCounter;

public:
	int index;
	int lastAnalyzed;

	Voxel();
	Voxel(double* da, double s, int i);
	Voxel(const Voxel & other);

	virtual ~Voxel();

	void miss();

	int getMissCounter();

	void resetMissCounter();

	bool isInside(double *da, double size);

	double *getPosition();


};

#endif /* VOXEL_H_ */
