/*
 * Parameters.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "Parameters.h"

Parameters::Parameters() {
	maxTriesWithoutSuccess = std::numeric_limits<int>::max();
	maxVoxels = 1000000;
	minDx = 0.0;
	from = 0;
	collectors = 1;
	maxTime = std::numeric_limits<double>::infinity();
	analyze = 10;
	split = 3000;
	surfaceSize = pow(100.0, 1.0/2.0);

	particleType = "Sphere";
	particleAttributes = "2";
}

Parameters::~Parameters() {
	// TODO Auto-generated destructor stub
}

