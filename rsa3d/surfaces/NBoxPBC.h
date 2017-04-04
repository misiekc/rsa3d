/*
 * NBoxPBC.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef SURFACES_NBOXPBC_H_
#define SURFACES_NBOXPBC_H_

#include "../Surface.h"
#include "../RND.h"

class NBoxPBC : public Surface {
private:
	static double * getTranslation(double *result, double s, double *p1, double *p2);

public:
	NBoxPBC(int dim, double s, double ndx, double vdx);
	virtual ~NBoxPBC();

	double getArea();
	double * getTranslation(double *result, double *p1, double *p2);
	void vector(double *v);
	double * getRandomPosition(double * result, RND rnd);
};

#endif /* SURFACES_NBOXPBC_H_ */
