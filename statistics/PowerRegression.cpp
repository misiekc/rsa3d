/*
 * LinearRegression.cpp
 *
 *  Created on: 22.06.2017
 *      Author: ciesla
 */

#include "PowerRegression.h"

#include <algorithm>
#include <cmath>

PowerRegression::PowerRegression() {
	this->s2a = 0.0;
	this->a = 0.0;
	this->b = 0.0;
}

PowerRegression::~PowerRegression() {
	for(DataElement *d: this->data)
		delete d;
}


/**
 * Add data point
 * @param x
 * @param y
 */
void PowerRegression::addXY(double x, double y){
	DataElement* d = new DataElement();
	d->x = x;
	d->y = y;
	this->data.push_back(d);
}

/**
 * calculates fit values
 */
void PowerRegression::calculate(int from, int to){
	double slnx = 0.0;
	double slny = 0.0;
	double slnxlny = 0.0;
	double slnx2 = 0.0;
	for (int i=from; i<to; i++){
		DataElement *d = this->data[i];
		double lnx = log(d->x);
		double lny = log(d->y);
		slnxlny += lnx*lny;
		slnx += lnx;
		slny += lny;
		slnx2 += lnx*lnx;
	}
	this->a = ((to-from)*slnxlny - slnx*slny) / ((to-from)*slnx2 - slnx*slnx);
	this->b = (slny - this->a*slnx)/(to-from);

	double s2 = 0.0;
	for (int i=from; i<to; i++){
		DataElement *d = this->data[i];
		s2 += pow(log(d->y) - this->b - this->a*log(d->x), 2.0);
	}
	this->s2a = ((to-from)/((to-from)*slnx2 - slnx*slnx))*((s2)/(to-from-2));
}

/**
 * calculates fit values
 */
void PowerRegression::calculate(){
	this->calculate(0, this->data.size());
}

/**
 * @return parameter A from y = Bx^A
 */
double PowerRegression::getA(){
	return this->a;
}

/**
 * @return standard deviation squared for parameter A from y = Ax + B
 */
double PowerRegression::getSA(){
	return sqrt(this->s2a);
}

/**
 * @return parameter B from y = Bx^A
 */
double PowerRegression::getB(){
	return exp(this->b);
}

int PowerRegression::size(){
	return this->data.size();
}

/**
 * Tester
 * @param args
 */
/*
void main(int argc, char **argv){
	double x, y;
	LinearRegression lr = new LinearRegression();
	for (int i=0; i<100; i++){
		x = i;
		y = (2+0.5*(Math.random()-0.5))*x + 5  + (Math.random()-0.5);
		System.out.println(x + " " + y);
		lr.addXY(x, y);
	}
	lr.calculate();
	System.out.println(
			String.valueOf(lr.getA()) + "\t" + //$NON-NLS-1$
			String.valueOf(lr.getSA()) + "\t" + //$NON-NLS-1$
			String.valueOf(lr.getB()) + "\t" + //$NON-NLS-1$
			String.valueOf(lr.getSB()) + "\t" + //$NON-NLS-1$
			String.valueOf(lr.getR())
			);
	}
	delete lr;
}
*/
