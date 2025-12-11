/*
 * Plot.cpp
 *
 *  Created on: 26.06.2017
 *      Author: ciesla
 */

#include "Plot.h"
#include <cmath>
#include <algorithm>

Plot::Plot(double min, double max, size_t bins) {

	this->yValues = new double[bins];
	this->y2Values = new double[bins];
	this->yCounter = new size_t[bins];

	for(size_t i=0; i<bins; i++){
		this->yValues[i] = 0.0;
		this->y2Values[i] = 0.0;
		this->yCounter[i] = 0;

	}
	this->bins = bins;

	this->min = min;
	this->max = max;
	this->step = (this->max - this->min) / static_cast<double>(bins);
}

Plot::~Plot() {
	delete[] this->yValues;
	delete[] this->y2Values;
	delete[] this->yCounter;
}

double Plot::getMax() const {
	return this->max;
}

size_t Plot::size() const {
	return this->bins;
}

size_t Plot::getIndex(double x) {
    if (x < this->min) {
        return 0;
    }
	return static_cast<size_t>((x - this->min) / this->step);
}

/**
 * adds point to plot (histogram)
 *
 * @param x
 */
void Plot::add(double x) {
	size_t i = this->getIndex(x);
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
	size_t i = this->getIndex(x);
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
	size_t index0 = this->getIndex(x0);
	size_t index1 = this->getIndex(x1);
//	if (index0 == index1) {
//		this.add(x0, y);
//	} else {
		for (size_t i = index0; i < index1; i++) {
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
	for (size_t i = 0; i < this->bins; i++) {
		x = this->min + static_cast<double>(i) * this->step + this->step / 2;
		if (this->yCounter[i]!=0)
			y = this->yValues[i] / static_cast<double>(this->yCounter[i]);
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
	for (size_t i = 0; i < this->bins; i++) {
		x = this->min + static_cast<double>(i) * this->step + this->step / 2;
		if (this->yCounter[i]!=0){
			y = this->yValues[i] / static_cast<double>(this->yCounter[i]);
			z = sqrt(this->y2Values[i] / static_cast<double>(this->yCounter[i]) - y*y);
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
	for (size_t i = 0; i < this->bins; i++) {
		x = this->min + static_cast<double>(i) * this->step + this->step / 2;
		y = static_cast<double>(this->yCounter[i]);
		points[i][0] = x;
		points[i][1] = y;
	}
	return points;
}

size_t Plot::getTotalNumberOfHistogramPoints(){
	size_t res = 0;
	for (size_t i = 0; i < this->bins; i++) {
		res += this->yCounter[i];
	}
	return res;
}
