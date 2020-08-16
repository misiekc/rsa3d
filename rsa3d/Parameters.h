/*
 * Parameters.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <limits>
#include <string>
#include <cmath>

#include "utils/Utils.h"
#include "utils/OMPMacros.h"

class Parameters {
private:
    void validateData();

public:
	std::size_t maxTriesWithoutSuccess = std::numeric_limits<std::size_t>::max();
	std::size_t maxVoxels = 40000000;
	double requestedAngularVoxelSize = 2*M_PI;
	double minDx = 0.0;
	std::size_t from = 0;
	std::size_t collectors = 1;
    double maxTime = std::numeric_limits<double>::infinity();
    double maxDensity = 1.0;
	std::size_t split = 5000;
	unsigned short surfaceDimension = RSA_SPATIAL_DIMENSION;
	double surfaceSize = pow(100000.0, 1.0/surfaceDimension);
	bool storePackings = true;
	std::string seedOrigin = "0";
	bool appendToDat = false;

	bool modifiedRSA = false;
	double thresholdDistance = 0.0;

	std::string boundaryConditions = "periodic";
	std::string particleType = "Sphere";
	std::string particleAttributes;

	std::size_t generatorProcesses = 1;
	int ompThreads = _OMP_MAXTHREADS;
	bool timestamp = false;

	bool coverageByNumber = false;

    bool operator==(const Parameters &rhs) const;
    bool operator!=(const Parameters &rhs) const;

    Parameters();
	explicit Parameters(std::istream &stream);
    explicit Parameters(const std::string &fileName);

    std::string getPackingSignature() const;
    double sufraceVolume() const;
};

#endif /* PARAMETERS_H_ */

