/*
 * ASFRegression.h
 *
 *  Created on: 30.06.2017
 *      Author: ciesla
 */

#ifndef ASFREGRESSION_H_
#define ASFREGRESSION_H_

#include <vector>

class ASFRegression {

	class DataElement{
	public:
		double x, y;
	};

private:
	std::vector<DataElement *> data;

	double s, sx, sx2, sx3, sx4, sxy, sx2y, delta;
	double c1, c2, sigmay, sigmac1, sigmac2;

public:
	ASFRegression();
	virtual ~ASFRegression();

	void addXY(double x, double y);
	void calculate(unsigned int from, unsigned int to);
	void calculate();
	void removeDistantPoints(double multiplier);
	double getC1();
	double getSC1();
	double getC2();
	double getSC2();
	int size();
};

#endif /* ASFREGRESSION_H_ */
