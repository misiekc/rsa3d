/*
 * Parameters.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "Parameters.h"
#include "Utils.h"
#include <iostream>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif

Parameters::Parameters() {
//	dimension = 2;
	maxTriesWithoutSuccess = std::numeric_limits<int>::max();
	maxVoxels = 40000000;
	minDx = 0.0;
	from = 0;
	collectors = 1;
	maxTime = std::numeric_limits<double>::infinity();
	split = 500;
	surfaceSize = pow(1000.0, 1.0/RSA_SPATIAL_DIMENSION);
	storePackings = true;

	modifiedRSA = false;
	thresholdDistance = 0.0;

	boundaryConditions = "periodic";
	particleType = "Sphere";
	particleAttributes = "2";
	generatorProcesses = 1;
#ifdef _OPENMP
	ompThreads =  omp_get_max_threads();
#else
	ompThreads = 1;
#endif

}

Parameters::Parameters(const std::string& sFile) : Parameters(){
	std::ifstream file(sFile, std::ios::in);
	if (!file) {
        std::cerr << "Cannot open configuration file: " + sFile << std::endl;
        exit(1);
    }

	std::string sLine;
	while( std::getline(file, sLine) ){
		if (sLine.length()==0 || sLine[0]=='#')
			continue;
		size_t pos;
		if( (pos=sLine.find("="))==std::string::npos)
			continue;
		std::string sKey = sLine.substr(0, pos);
		std::string sValue = sLine.substr(pos+1, sLine.length()-pos-1);
		trim(sKey);
		trim(sValue);

		if (sKey.compare("maxTriesWithoutSuccess")==0) 	this->maxTriesWithoutSuccess = std::stoi(sValue);
		else if (sKey.compare("maxVoxels")==0)					this->maxVoxels = std::stoi(sValue);
		else if (sKey.compare("requestedAngularVoxelSize")==0)	this->requestedAngularVoxelSize = std::stod(sValue);
		else if (sKey.compare("minDx")==0)						this->minDx = std::stod(sValue);
		else if (sKey.compare("maxTime")==0) 					this->maxTime = std::stod(sValue);
		else if (sKey.compare("split")==0) 						this->split = std::stoi(sValue);
		else if (sKey.compare("surfaceVolume")==0) 				this->surfaceSize = pow(std::stod(sValue), 1.0/RSA_SPATIAL_DIMENSION);
		else if (sKey.compare("storePackings")==0)	 			this->storePackings = (sValue.compare("false")==0)?false:true;
		else if (sKey.compare("modifiedRSA")==0)		 		this->modifiedRSA = (sValue.compare("false")==0)?false:true;
		else if (sKey.compare("thresholdDistance")==0) 			this->thresholdDistance = std::stod(sValue);
		else if (sKey.compare("boundaryConditions")==0) 		this->boundaryConditions = sValue;
		else if (sKey.compare("particleType")==0) 				this->particleType = sValue;
		else if (sKey.compare("particleAttributes")==0)			this->particleAttributes = sValue;

		else if (sKey.compare("from")==0) 						this->from = std::stoi(sValue);
		else if (sKey.compare("collectors")==0) 				this->collectors = std::stoi(sValue);
		else if (sKey.compare("generatorProcesses")==0) 		this->generatorProcesses = std::stoi(sValue);
		else if (sKey.compare("ompThreads")==0) 				this->ompThreads = std::stoi(sValue);
	}
	file.close();

	#ifdef _OPENMP
	omp_set_num_threads(this->ompThreads);
	#endif

}

Parameters::~Parameters() {
	// TODO Auto-generated destructor stub
}

