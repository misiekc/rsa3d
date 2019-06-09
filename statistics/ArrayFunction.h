/*
 * ArrayFunction.h
 *
 *  Created on: Jun 7, 2019
 *      Author: ciesla
 */

#ifndef STATISTICS_ARRAYFUNCTION_H_
#define STATISTICS_ARRAYFUNCTION_H_

class ArrayFunction {

private:
	double **data;
	unsigned int length;
public:
	ArrayFunction(double *data[2], unsigned int length);
	virtual ~ArrayFunction();

	double get(double x);
};

#endif /* STATISTICS_ARRAYFUNCTION_H_ */
