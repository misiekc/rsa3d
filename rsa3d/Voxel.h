/*
 * Voxel.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef VOXEL_H_
#define VOXEL_H_

class Voxel {
private:
	int dimension;
	double* center;
	int index;
	short missCounter;
	int lastAnalyzed;

public:
	Voxel(int dim, double* da, double s, int i);

	virtual ~Voxel();

	void miss();

	int getMissCounter();

	void resetMissCounter();

	double* getPosition();

};

#endif /* VOXEL_H_ */
