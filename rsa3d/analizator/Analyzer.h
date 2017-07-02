/*
 * Analyzer.h
 *
 *  Created on: 25.06.2017
 *      Author: ciesla
 */

#ifndef ANALIZATOR_ANALYZER_H_
#define ANALIZATOR_ANALYZER_H_

#include "../Shape.h"

#include <Plot.h>
#include <LogPlot.h>

#include <string>
#include <vector>

class Analyzer {
public:
	Analyzer(Parameters *params);
	virtual ~Analyzer();

	void analyzePackingsInDirectory(char *sdir, double mintime, double particleSize);

private:
	Parameters *params;

	std::vector<Shape *> * fromFile(std::string filename);
	void analyzePacking(std::vector<Shape *> *packing, LogPlot *nvt, Plot *asf, Plot *corr, double surfaceFactor);
	double * printNvT(LogPlot &nvt, std::string filename, double surfaceFactor, double *res);
	double * printASF(Plot &asf, std::string filename, int counter, double packingFraction, double *res);
	void printCorr(Plot &corr, std::string filename, int counter, double particleSize, double packingFraction);
};

#endif /* ANALIZATOR_ANALYZER_H_ */
