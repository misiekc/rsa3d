/*
 * ASFRegression.cpp
 *
 *  Created on: 30.06.2017
 *      Author: ciesla
 */

#include "ASFRegression.h"
#include <cmath>
#include <algorithm>

ASFRegression::ASFRegression() {
	this->s = 0.0;
	this->sx = 0.0;
	this->sx2 = 0.0;
	this->sx3 = 0.0;
	this->sx4 = 0.0;
	this->sxy = 0.0;
	this->sx2y = 0.0;
	this->c1 = 0.0;
	this->c2 = 0.0;
	this->sigmac1 = 0.0;
	this->sigmac2 = 0.0;
	this->sigmay = 0.0;
	this->delta = 0.0;
}

ASFRegression::~ASFRegression() {
	for(DataElement *d: this->data)
		delete d;
}

/**
 * Add data point
 * @param x
 * @param y
 */
void ASFRegression::addXY(double x, double y){
	DataElement* d = new DataElement();
	d->x = x;
	d->y = y;
	this->data.push_back(d);
}

/**
 * calculates fit values
 */
void ASFRegression::calculate(unsigned int from, unsigned int to){
		this->s = 0.0;
		this->sx = 0.0;
		this->sx2 = 0.0;
		this->sx3 = 0.0;
		this->sx4 = 0.0;
		this->sxy = 0.0;
		this->sx2y = 0.0;
		to = (to < this->data.size()) ? to : this->data.size();
		for (unsigned int i=from; i<to; i++){
			DataElement *d = this->data[i];
			this->sx += d->x;
			this->sx2 += d->x*d->x;
			this->sx3 += d->x*d->x*d->x;
			this->sx4 += d->x*d->x*d->x*d->x;

			this->sxy += d->x*d->y;
			this->sx2y += d->x*d->x*d->y;
		}

		this->delta = this->sx3*this->sx3 - this->sx2*this->sx4;

		this->c1 = -(-this->sx2*this->sx3 + this->sx*this->sx4 - this->sx4*this->sxy + this->sx3*this->sx2y) / this->delta;

		this->c2 = -(-this->sx2*this->sx2 + this->sx*this->sx3 - this->sx3*this->sxy + this->sx2*this->sx2y) / this->delta;


		this->sigmay = 0.0;
		for(unsigned int i=from; i<to; i++){
			DataElement *d = this->data[i];
			this->sigmay += pow( d->y - (1 - this->c1*d->x + this->c2*d->x*d->x), 2.0);
		}
		this->sigmay /= (to-from-2);

		this->sigmac1 = sqrt(this->sigmay)* fabs((this->sx4*this->sx - this->sx3*this->sx2) / this->delta);
		this->sigmac2 = sqrt(this->sigmay)* fabs((this->sx3*this->sx - this->sx2*this->sx2) / this->delta);
	}

/**
 * calculates fit values
 */
void ASFRegression::calculate(){
	this->calculate(0, this->data.size());
}

/**
 * @param multiplier
 */
void ASFRegression::removeDistantPoints(double multiplier){
	double sigma = this->sigmay/this->data.size();
	for (DataElement *d : this->data) {
		if (pow(d->y - (1 - this->c1*d->x + this->c2*d->x*d->x), 2.0) > multiplier*sigma){
			this->data.erase(std::remove(this->data.begin(), this->data.end(), d), this->data.end());
		}
	}
	this->calculate();
}


/**
 * @return parameter A from y = Ax + B
 */
double ASFRegression::getC1(){
	return this->c1;
}

/**
 * @return standard deviation squared for parameter A from y = Ax + B
 */
double ASFRegression::getSC1(){
	return this->sigmac1;
}

/**
 * @return parameter B from y = Ax + B
 */
double ASFRegression::getC2(){
	return this->c2;
}

/**
 * @return standard deviation sqared for parameter B from y = Ax + B
 */
double ASFRegression::getSC2(){
	return this->sigmac2;
}

int ASFRegression::size(){
	return this->data.size();
}
