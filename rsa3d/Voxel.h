/*
 * Voxel.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef VOXEL_H_
#define VOXEL_H_

#include "Positioned.h"

class Voxel : public Positioned{

private:
	short missCounter;

public:
	int index;
	int lastAnalyzed;

	Voxel(unsigned char dim);
	Voxel(unsigned char dim, double* da, double s, int i);
	Voxel(const Voxel & other);

	virtual ~Voxel();

	void miss();

	int getMissCounter();

	void resetMissCounter();


};

#endif /* VOXEL_H_ */
