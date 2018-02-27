/*
 * Voxel.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef VOXEL_H_
#define VOXEL_H_

#include <algorithm>

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
class Voxel{

private:
	double position[SPATIAL_DIMENSION];
	double orientation[ANGULAR_DIMENSION];
	short missCounter;

public:
	int index;
	int lastAnalyzed;

	Voxel();
	Voxel(double* pos, double *angle, int i);
	Voxel(const Voxel & other);

	virtual ~Voxel();

	void miss();

	int getMissCounter();

	void resetMissCounter();

	bool isInside(double *pos, double size);

	bool isInside(double *pos, double size, double* angle, double asize);

	double *getPosition();

	double *getOrientation();


};

#include "Voxel.tpp"

#endif /* VOXEL_H_ */
