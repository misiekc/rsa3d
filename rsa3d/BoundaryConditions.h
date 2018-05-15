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

	virtual double distance2(const double *p1, const double *p2) = 0;
	/**
	 * @brief Returns translation that should be applied to @a p2 to move him to the "proximity" of @a p1
	 */
	virtual double * getTranslation(double *result, const double* p1, const double* p2) = 0;

};

#endif /* BOUNDARYCONDITIONS_H_ */
