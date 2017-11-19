/*
 * Voxel.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef VOXEL_H_
#define VOXEL_H_

#include <algorithm>

template <ushort DIMENSION>
class Voxel{

private:
	double position[DIMENSION];
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

#include "Voxel.tpp"

#endif /* VOXEL_H_ */
