/*
 * LogPlot.h
 *
 *  Created on: 27.06.2017
 *      Author: ciesla
 */

#ifndef LOGPLOT_H_
#define LOGPLOT_H_

#include "Plot.h"

class LogPlot: public Plot {
private:
	double factor;

protected:
	int getIndex(double x) override;

public:
	LogPlot(double min, double max, int bins);
	LogPlot(double max, int bins);
	virtual ~LogPlot();

	virtual double** getAsPoints(double **points) override;
	virtual double** getAsPointsWithErrors(double **points) override;
};

#endif /* LOGPLOT_H_ */
