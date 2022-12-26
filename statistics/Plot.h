/*
 * Plot.h
 *
 *  Created on: 26.06.2017
 *      Author: ciesla
 */

#ifndef PLOT_H_
#define PLOT_H_

#include <cstdlib>

class Plot {
protected:
	double* yValues;
	double* y2Values;
	size_t* yCounter;
	size_t bins;
	double min, max, step;

	virtual size_t getIndex(double x);

public:
	Plot(double min, double max, size_t bins);
	Plot(const Plot &other) = delete;
	Plot &operator=(const Plot &other) = delete;

	virtual ~Plot();
	[[nodiscard]] double getMax() const;
	[[nodiscard]] size_t size() const;
	void add(double x);
	void add(double x, double y);
	void addBetween(double x0, double x1, double y);
	virtual double** getAsPoints(double **points);
	virtual double** getAsPointsWithErrors(double **points);
	virtual double** getAsHistogramPoints(double **points);
	size_t getTotalNumberOfHistogramPoints();
};

#endif /* PLOT_H_ */
