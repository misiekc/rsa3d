/*
 * Analyzer.h
 *
 *  Created on: 25.06.2017
 *      Author: ciesla
 */

#ifndef ANALIZATOR_ANALYZER_H_
#define ANALIZATOR_ANALYZER_H_

#include "../Shape.h"
#include "../Parameters.h"
#include "../shapes/Cuboid.h"

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

	std::vector<Shape<RSA_DIMENSION> *> * fromFile(std::string filename);
	void analyzePacking(std::vector<Shape<RSA_DIMENSION> *> *packing, LogPlot *nvt, Plot *asf, Plot *corr, double surfaceFactor);
	void analyzeOrder(std::vector<Shape<RSA_DIMENSION> *> *packing, Plot **order);
	void calculateOrderParameters(double *result, Cuboid *c1, Cuboid *c2);
	double * printNvT(LogPlot &nvt, std::string filename, double *fixedA, double surfaceFactor, double *res);
	double * printASF(Plot &asf, std::string filename, int counter, double packingFraction, double *res);
	void printCorrelations(Plot &corr, std::string filename, int counter, double particleSize, double packingFraction);
	void printOrder(Plot **order, std::string filename);
};

#endif /* ANALIZATOR_ANALYZER_H_ */
