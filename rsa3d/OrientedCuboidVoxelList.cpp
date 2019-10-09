/*
 * OrientedCuboidVoxelList.cpp
 *
 *  Created on: 30.09.2019
 *      Author: ciesla
 */

#include <iostream>
#include <cmath>
#include <algorithm>

#include "OrientedCuboidVoxelList.h"
#include "utils/OMPMacros.h"
#include "utils/Assertions.h"

#include "shape/shapes/OrientedCuboid.h"

OrientedCuboidVoxelList::OrientedCuboidVoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange, double requestedAngularVoxelSize) : VoxelList(dim, packingSpatialSize, requestedSpatialVoxelSize, shapeAngularRange, requestedAngularVoxelSize){
}

/**
 * returns true when the whole voxel is inside an exclusion area of any shape in shapes
 * To determine it the method tires to split voxel up to level of maxDepth
 */
bool OrientedCuboidVoxelList::isVoxelInsideExclusionZone(Voxel *v, double spatialSize, double angularSize,
										   std::vector<const RSAShape *> *shapes, RSABoundaryConditions *bc,
                                           unsigned short depth){

	size_t finalArrayLength = (size_t)round( pow(pow(2.0, depth), this->surfaceDimension+RSA_ANGULAR_DIMENSION) );
	Voxel **finalVoxels = new Voxel*[ finalArrayLength ];
	double ss = spatialSize;
	double as = angularSize;
	size_t length = 1;
	// dzielimy voxel na mniejsze (gdy depth > 0)
	size_t tmpArrayLength = (size_t)round( pow(2, this->surfaceDimension+RSA_ANGULAR_DIMENSION) );
	Voxel **tmpVoxels = new Voxel*[ tmpArrayLength ];
	finalVoxels[0] = v;
	size_t last = 0;
	for(unsigned short i=0; i<depth; i++){
		ss /= 2.0;
		as /= 2.0;
		for(size_t j=0; j<length; j++){
			this->splitVoxel(finalVoxels[j], ss, as, tmpVoxels);
			finalVoxels[j] = tmpVoxels[0];
			for(size_t k=1; k<tmpArrayLength; k++){
				last++;
				finalVoxels[last] = tmpVoxels[k];
			}
		}
		length *= tmpArrayLength;
	}
	Validate(last==length-1);
	// sprawdzamy voxele z tablicy finalVoxels.
	bool allInside = true;
	for(size_t i=0; i<finalArrayLength; i++){
		// szukamy ksztaltu, ktorego strefa wykluczjaca zawiera woksel
		bool isInside = false; // będziemy szukać wartosci true
		isInside = OrientedCuboid<RSA_SPATIAL_DIMENSION>::voxelInside(bc, finalVoxels[i]->getPosition(), ss, shapes);
/*
		for(const RSAShape *s : *shapes){
			// jesli choc jeden z mniejszych jest false zwracamy false
			isInside = s->voxelInside(bc, finalVoxels[i]->getPosition(), finalVoxels[i]->getOrientation(), ss, as);
			if (isInside)
				break; // znaleziony, przechodzimy do nastepnego woksela
		}
*/
		if (!isInside){
			allInside = false;
			break; // zaden ksztalt nie zawiera tego woksela. zwracamy false
		}
	}
	if (depth>0){
		for(size_t i=0; i<finalArrayLength; i++){
			delete finalVoxels[i];
		}
	}
	delete[] finalVoxels;
	delete[] tmpVoxels;
	return allInside;
}
