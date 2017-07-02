/*
 * Plot.h
 *
 *  Created on: 26.06.2017
 *      Author: ciesla
 */

#ifndef PLOT_H_
#define PLOT_H_

class Plot {
protected:
	double* yValues;
	double* y2Values;
	int* yCounter;
	int bins;
	double min, max, step;

	virtual int getIndex(double x);

public:
	Plot(double min, double max, int bins);
	virtual ~Plot();
	double getMax();
	int size();
	void add(double x, double y);
	void addBetween(double x0, double x1, double y);
	virtual double** getAsPoints(double **points);
	virtual double** getAsPointsWithErrors(double **points);

	double** getAsHistogramPoints(double **points);
};

#endif /* PLOT_H_ */
