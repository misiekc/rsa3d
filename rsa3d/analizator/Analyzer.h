/*
 * Analyzer.h
 *
 *  Created on: 25.06.2017
 *      Author: ciesla
 */

#ifndef ANALIZATOR_ANALYZER_H_
#define ANALIZATOR_ANALYZER_H_

#include "../shape/Shape.h"
#include "../Parameters.h"
#include "../shape/shapes/cuboid/Cuboid.h"

#include "../../statistics/Plot.h"
#include "../../statistics/LogPlot.h"
#include "../Packing.h"

#include <string>
#include <vector>

class Analyzer {
public:
	explicit Analyzer(Parameters *params) : params(params) {}
	virtual ~Analyzer() = default;

	void analyzePackingsInDirectory(char *sdir, double mintime, double particleSize);

private:
	Parameters *params;

	void analyzePacking(const Packing &packing, LogPlot *nvt, Plot *asf, Plot *corr, double surfaceFactor);
	void analyzeOrder(const Packing &packing, std::vector<Plot*> *order);
	double * printNvT(LogPlot &nvt, std::string filename, double *fixedA, double surfaceFactor, double *res);
	double * printASF(Plot &asf, std::string filename, int counter, double packingFraction, double *res);
	void printCorrelations(Plot &corr, std::string filename, int counter, double particleSize, double packingFraction);
	void printOrder(const std::vector<Plot*> &order, const std::string &filename) const;
//	double getPetiodicDistance(const RSAShape *shape1, const RSAShape *shape2) const;
	std::vector<Plot *> getFilledOrderVector() const;
	bool isOrderCalculable(const RSAShape *shape) const;
};

#endif /* ANALIZATOR_ANALYZER_H_ */
