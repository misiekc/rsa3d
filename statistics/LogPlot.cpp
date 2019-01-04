/*
 * LogPlot.cpp
 *
 *  Created on: 27.06.2017
 *      Author: ciesla
 */

#include "LogPlot.h"
#include <cmath>

LogPlot::LogPlot(double min, double max, int bins) : Plot(min, max, bins){
	this->factor = pow(this->max/this->min, 1.0/bins);
}

LogPlot::LogPlot(double max, int bins) : LogPlot(1.0, max, bins){
}

LogPlot::~LogPlot() {
}


/*
	public LogPlot(Histogram h, boolean decount, boolean normalize){
		super(h.getBins()/100);
		this.min = 1;
		this.max = h.getMax();
		this.factor = Math.pow(this.max/this.min, 1.0/(h.getBins()/100));
		Point2D[] points = h.getAsPoints(decount, normalize);
		for(int i=0; i<points.length; i++){
			this.add(points[i].getX(), points[i].getY());
		}
	}
*/

int LogPlot::getIndex(double d){
	d /= this->min;
	if (d>1.0)
		return (int)(log(d)/log(this->factor));
	return 0;
}

/**
 * @return plot as point array
 */
double ** LogPlot::getAsPoints(double **points){
	double x, y;
	int i;
	for(i=0; i<this->bins; i++){
		x = this->min*pow(this->factor, i+0.5);
		if (this->yCounter[i]!=0)
			y = this->yValues[i] / this->yCounter[i];
		else
			y = 0;
		points[i][0] = x;
		points[i][1] = y;
	}
	return points;
}

	/**
	 * @return plot as point array
	 */
double** LogPlot::getAsPointsWithErrors(double **points){
	double x, y, z;
	int i;
	for(i=0; i<this->bins; i++){
		x = this->min*pow(this->factor, i+0.5);
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
