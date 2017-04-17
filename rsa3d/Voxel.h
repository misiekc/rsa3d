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

	Voxel(int dim);
	Voxel(int dim, double* da, double s, int i);

	virtual ~Voxel();

	void miss();

	int getMissCounter();

	void resetMissCounter();


};

#endif /* VOXEL_H_ */
