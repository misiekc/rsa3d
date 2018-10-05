/*
 * Parameters.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "Parameters.h"
#include "Utils.h"
#include "Config.h"
#include <iostream>
#include <fstream>


Parameters::Parameters(std::istream &stream) {
	auto config = Config::parse(stream);
	for (const auto &key : config.getKeys()){
		if (key == "maxTriesWithoutSuccess") 	        this->maxTriesWithoutSuccess = config.getUnsignedLong(key);
		else if (key == "maxVoxels")					this->maxVoxels = config.getUnsignedLong(key);
		else if (key == "requestedAngularVoxelSize")	this->requestedAngularVoxelSize = config.getDouble(key);
		else if (key == "minDx")						this->minDx = config.getDouble(key);
		else if (key == "maxTime") 					    this->maxTime = config.getDouble(key);
		else if (key == "split") 						this->split = config.getUnsignedLong(key);
		else if (key == "surfaceVolume") 				this->surfaceSize = pow(config.getDouble(key), 1.0/RSA_SPATIAL_DIMENSION);
		else if (key == "storePackings")	 			this->storePackings = config.getString(key) != "false";
		else if (key == "modifiedRSA")		 		    this->modifiedRSA = config.getString(key) != "false";
		else if (key == "thresholdDistance") 			this->thresholdDistance = config.getDouble(key);
		else if (key == "boundaryConditions") 		    this->boundaryConditions = config.getString(key);
		else if (key == "particleType") 				this->particleType = config.getString(key);
		else if (key == "particleAttributes")			this->particleAttributes = config.getString(key);

		else if (key == "from") 						this->from = config.getUnsignedLong(key);
		else if (key == "collectors") 				    this->collectors = config.getUnsignedLong(key);
		else if (key == "generatorProcesses") 		    this->generatorProcesses = config.getUnsignedLong(key);
		else if (key == "ompThreads") 				    this->ompThreads = config.getInt(key);
	}

	this->validateData();

	#ifdef _OPENMP
	omp_set_num_threads(this->ompThreads);
	#endif
}

Parameters::Parameters(const std::string &fileName) {
    std::ifstream input(fileName, std::ios::in);
    if (!input)
        throw std::runtime_error("Cannot open configuration file: " + fileName);

    (*this) = Parameters{input};    // Yes, it is possible ;)
}

void Parameters::validateData() {
    Ensure(maxTriesWithoutSuccess > 0);
    Ensure(maxVoxels > 0);
    Ensure(requestedAngularVoxelSize > 0);
    Ensure(minDx >= 0.0);
    Ensure(from >= 0);
    Ensure(collectors > 0);
    Ensure(maxTime > 0);
    Ensure(split > 0);
    Ensure(!std::isnan(surfaceSize) && surfaceSize > 0);
    Ensure(thresholdDistance >= 0.0);
    Ensure(generatorProcesses > 0);
    Ensure(ompThreads > 0);
    Ensure(!particleType.empty());
}

Parameters::Parameters() {
    this->validateData();       // Make sure default construction is valid
}
