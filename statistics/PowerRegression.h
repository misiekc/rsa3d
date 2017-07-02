/*
 * PowerRegression.h
 *
 *  Created on: 22.06.2017
 *      Author: ciesla
 */

#ifndef POWERREGRESSION_H_
#define POWERREGRESSION_H_

#include <vector>

class PowerRegression {
	class DataElement{
	public:
		double x, y;
	};

private:

	std::vector<DataElement *> data;

	double a, b, s2a;

public:
	PowerRegression();
	virtual ~PowerRegression();

	void addXY(double x, double y);
	void calculate(int from, int to);
	void calculate();
	void removeDistantPoints(double multiplier);
	double getA();
	double getSA();
	double getB();
	int size();
};

#endif /* ANALIZATOR_LINEARREGRESSION_H_ */
