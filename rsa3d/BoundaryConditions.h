/*
 * BoundaryConditions.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef BOUNDARYCONDITIONS_H_
#define BOUNDARYCONDITIONS_H_

class BoundaryConditions {
public:
	BoundaryConditions();
	virtual ~BoundaryConditions();

	virtual double distance2(double *p1, double *p2);
	virtual double* getTranslation(double* p1, double* p2);

};

#endif /* BOUNDARYCONDITIONS_H_ */
