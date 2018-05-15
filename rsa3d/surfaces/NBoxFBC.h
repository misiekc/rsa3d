/*
 * NBoxFBC.h
 *
 *  Created on: 13.07.2017
 *      Author: ciesla
 */

#ifndef SURFACES_NBOXFBC_H_
#define SURFACES_NBOXFBC_H_

#include "../Surface.h"
#include "../RND.h"

class NBoxFBC: public Surface {
private:
	static double * getTranslation(double *result, int dim, double s, const double *p1, const double *p2);

public:
	NBoxFBC(int dim, double s, double ndx, double vdx);
	virtual ~NBoxFBC();

	double getArea();
	double * getTranslation(double *result, const double *p1, const double *p2);
	void vector(double *v);
};

#endif /* SURFACES_NBOXFBC_H_ */
