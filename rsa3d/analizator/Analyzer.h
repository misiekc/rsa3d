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
#include "../PackingGenerator.h"
#include "../NeighbourGrid.h"
#include "../utils/Quantity.h"


#include <string>
#include <vector>

class Analyzer {
public:
	explicit Analyzer(const Parameters *params) : params(params) {}
	virtual ~Analyzer() = default;

	void analyzePackingsInDirectory(const std::string &dirName, double mintime, double particleSize,
                                    double correlationsRange);

private:
	const Parameters *params;

	struct Result {
		std::string dir;
		Quantity theta;
		Quantity d;
		Quantity thetaInf;
		Quantity C1;
		Quantity C2;
        Quantity A;

		void print(std::ostream &out) const;
    };

	std::vector<std::vector<double>> orderParameterRange;

	/* returns packing fraction */
	double analyzePacking(const Packing &packing, LogPlot *nvt, Plot *asf, double surfaceFactor);
	void analyzeCorrelations(const Packing &packing, const NeighbourGrid<const RSAShape> &ng, Plot *corr);
	const RSAShape *getNearestNeighbour(const RSAShape *s, const NeighbourGrid<const RSAShape> &ng);
	void analyzeNearestNeighbours(const Packing &packing, const NeighbourGrid<const RSAShape> &ng, Plot *minDist);
	void analyzeOrder(const Packing &packing, const NeighbourGrid<const RSAShape> &ng, std::vector<Plot*> *order);
	void printKinetics(LogPlot &nvt, std::string filename, double *fixedA, double surfaceFactor, double dTime, Result *res);
	void printASF(Plot &asf, std::string filename, int counter, double packingFraction, Result *res);
	void printHistogram(Plot& plot, std::string filename);
	void printCorrelations(Plot &corr, std::string filename);
	void fillOrderParameterRange(const RSAShape *s);
	double applyOrderNormalisation(double d, std::vector<double> &norm);
	void printOrder(const std::vector<Plot*> &order, const std::string &filename) const;
//	double getPetiodicDistance(const RSAShape *shape1, const RSAShape *shape2) const;
	std::vector<Plot *> getFilledOrderVector(double range) const;
	bool isOrderCalculable(const RSAShape *shape) const;

    void findMinMaxTimes(double *minmaxTimes, const std::vector<std::string> &packingPaths);
};

#endif /* ANALIZATOR_ANALYZER_H_ */
