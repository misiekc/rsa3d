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
	RSAOrientation requestedAngularVoxelSize;
    RSAOrientation angularVoxelRange;
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

	bool modifiedRSA = false;
	double thresholdDistance = 0.0;

	std::string boundaryConditions = "periodic";
	std::string surfaceFunction = "";
	std::string particleType = "Sphere";
	std::string particleAttributes;

	std::size_t generatorProcesses = 1;
	int ompThreads = _OMP_MAXTHREADS;
	unsigned short maxDeepLevel = 0;
	bool timestamp = false;

	bool coverageByNumber = false;

    bool operator==(const Parameters &rhs) const;
    bool operator!=(const Parameters &rhs) const;

    Parameters();
	explicit Parameters(std::istream &stream);
    explicit Parameters(const std::string &fileName);

    [[nodiscard]] std::string getPackingSignature() const;

    /**
     * @brief Returns surface volume assuming Parameters::surfaceSize linear size and
     * Parameters::packingFractionSurfaceDimension() dimension
     * @return
     */
    [[nodiscard]] double sufraceVolume() const;

    /**
     * @brief Return the dimension of surfaces with respect to which the packing fractio nshould be calculated.
     * @details For example, for curved surfaces, although the actual surface is n-dimentional, the packing fraction
     * is calculated a a projection on n-1 first dimensions.
     */
    [[nodiscard]] std::size_t packingFractionSurfaceDimension() const;
};

#endif /* PARAMETERS_H_ */

