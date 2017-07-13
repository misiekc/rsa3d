/*
 * Parameters.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <climits>
#include <limits>
#include <string>
#include <cmath>

class Parameters {

public:
	int dimension = 2;
	int maxTriesWithoutSuccess = std::numeric_limits<int>::max();
	int maxVoxels = 100000000;
	double minDx = 0.0;
	int from = 0;
	int collectors = 1;
	double maxTime = std::numeric_limits<double>::infinity();
	int analyze = 10;
	int split = 3000;
	double surfaceSize = pow(100000.0, 1.0/2.0);

	std::string boundaryConditions;
	std::string particleType;
	std::string particleAttributes;


	Parameters();
	Parameters(const std::string & sFile);
	virtual ~Parameters();
};

#endif /* PARAMETERS_H_ */

