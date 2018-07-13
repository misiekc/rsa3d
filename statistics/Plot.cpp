/*
 * Plot.cpp
 *
 *  Created on: 26.06.2017
 *      Author: ciesla
 */

#include "Plot.h"
#include <cmath>
#include <algorithm>

Plot::Plot(double min, double max, int bins) {

	this->yValues = new double[bins];
	this->y2Values = new double[bins];
	this->yCounter = new unsigned long int[bins];

	for(int i=0; i<bins; i++){
		this->yValues[i] = 0.0;
		this->y2Values[i] = 0.0;
		this->yCounter[i] = 0;

	}
	this->bins = bins;

	this->min = min;
	this->max = max;
	this->step = (this->max - this->min) / bins;
}

Plot::~Plot() {
	delete[] this->yValues;
	delete[] this->y2Values;
	delete[] this->yCounter;
}

double Plot::getMax() {
	return this->max;
}

int Plot::size() {
	return this->bins;
}

int Plot::getIndex(double x) {
	return (int) ((x - this->min) / this->step);
}

/**
 * adds point to plot (histogram)
 *
 * @param x
 */
void Plot::add(double x) {
	int i = this->getIndex(x);
	if (i >= 0 && i < this->bins) {
		this->yValues[i] += 0.0;
		this->y2Values[i] += 0.0;
		this->yCounter[i]++;
	}
}

/**
 * adds point to plot. Values added for the same x will be averaged
 *
 * @param x
 * @param y
 */
void Plot::add(double x, double y) {
	int i = this->getIndex(x);
	if (i >= 0 && i < this->bins) {
		this->yValues[i] += y;
		this->y2Values[i] += y*y;
		this->yCounter[i]++;
	}
}

/**
 * adds point to plot. Values added for the same x will be averaged
 *
 * @param x
 * @param y
 */
void Plot::addBetween(double x0, double x1, double y) {
	int index0 = this->getIndex(x0);
	int index1 = this->getIndex(x1);
//	if (index0 == index1) {
//		this.add(x0, y);
//	} else {
		for (int i = std::max(0, index0); i < index1; i++) {
			if (i < this->bins){
				this->yValues[i] += y;
				this->y2Values[i] += y*y;
				this->yCounter[i]++;
			}
		}
//	}
}

/**
 * @return plot as point array
 */
double** Plot::getAsPoints(double** points) {
	double x, y;
	int i;
	for (i = 0; i < this->bins; i++) {
		x = this->min + i * this->step + this->step / 2;
		if (this->yCounter[i]!=0)
			y = this->yValues[i] / this->yCounter[i];
		else
			y = 0.0;
		points[i][0] = x;
		points[i][1] = y;
	}
	return points;
}

/**
 * @return plot as point array
 */
double** Plot::getAsPointsWithErrors(double **points) {
	double x, y, z;
	int i;
	for (i = 0; i < this->bins; i++) {
		x = this->min + i * this->step + this->step / 2;
		if (this->yCounter[i]!=0){
			y = this->yValues[i] / this->yCounter[i];
			z = sqrt(this->y2Values[i] / this->yCounter[i] - y*y);
		}else{
			y = 0.0;
			z = 0.0;
		}
		points[i][0] = x;
		points[i][1] = y;
		points[i][2] = z;
	}
	return points;
}

double** Plot::getAsHistogramPoints(double** points) {
	double x, y;
	int i;
	for (i = 0; i < this->bins; i++) {
		x = this->min + i * this->step + this->step / 2;
		y = this->yCounter[i];
		points[i][0] = x;
		points[i][1] = y;
	}
	return points;
}

unsigned long int Plot::getTotalNumberOfHistogramPoints(){
	unsigned long int res = 0;
	for (int i = 0; i < this->bins; i++) {
		res += this->yCounter[i];
	}
	return res;
}
