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
	size_t getIndex(double x) override;

public:
	LogPlot(double min, double max, int bins);
	LogPlot(double max, int bins);
	LogPlot(const LogPlot &other) = delete;
	LogPlot &operator=(const LogPlot &other) = delete;
	~LogPlot() override;

	double** getAsPoints(double **points) override;
	double** getAsPointsWithErrors(double **points) override;
    double** getAsHistogramPoints(double** points) override;

    };

#endif /* LOGPLOT_H_ */
