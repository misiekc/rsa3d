/*
 * DiscreteAngleVoxelList.cpp
 *
 *  Created on: 11.11.2022
 *      Author: Michal Ciesla
 */

#include <iostream>
#include <cmath>
#include <algorithm>

#include "DiscreteAngleVoxelList.h"
#include "utils/OMPMacros.h"
#include "utils/Assertions.h"


DiscreteAngleVoxelList::DiscreteAngleVoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, const std::vector<RSAOrientation> &orientations) : VoxelList(){
	Expects(dim > 0);
	Expects(dim <= RSA_SPATIAL_DIMENSION);
	Expects(packingSpatialSize > 0.0);
	Expects(requestedSpatialVoxelSize > 0.0);

	this->surfaceDimension = dim;
	this->spatialRange = packingSpatialSize;
	this->allowedOrientations = orientations;


	this->initialVoxelSize = VoxelList::findFloorSize(requestedSpatialVoxelSize);
	this->spatialVoxelSize = this->spatialRange;
	for (unsigned short int i=0; i<RSA_ANGULAR_DIMENSION; i++)
		this->angularVoxelSize[i] = 0;
	for (unsigned short i=0; i<RSA_ANGULAR_DIMENSION; i++)
		this->angularRange[i] = 2*M_PI;
	this->activeTopLevelVoxels = nullptr;
	this->voxelNeighbourGrid = nullptr;


	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);
	this->disabled = false;
	this->voxelsInitialized = false;

	RSAVector position;
	for(unsigned short i=0; i<RSA_SPATIAL_DIMENSION; i++)
		position[i] = 0.0;

	this->voxels = new Voxel*[this->allowedOrientations.size()];
	for (size_t i=0; i<this->allowedOrientations.size(); i++){
		this->voxels[i] = new Voxel(position, this->allowedOrientations[i]);
		this->voxels[i]->index = i;
	}
	this->length = this->allowedOrientations.size();
}

unsigned int DiscreteAngleVoxelList::initVoxels(RSABoundaryConditions *bc, NeighbourGrid<const RSAShape> *nl){
    // "Full spatial" variables correspond to the whole cube simulation box (for initial voxels and neighbour grid)
    // "Spatial" (without "full") variables correspond to the acutal part where voxels will be created (for ordinary
    // voxels)
	size_t fullSpatialGridLinearSize = this->findArraySize(this->spatialRange, this->initialVoxelSize);
	size_t fullSpatialGridSize = (size_t)(pow(fullSpatialGridLinearSize, this->surfaceDimension) + 0.5);
	std::array<std::size_t, RSA_SPATIAL_DIMENSION> spatialGridLinearSize = this->calculateSpatialGridLinearSize();
	Assert(RSA_SPATIAL_DIMENSION >= this->surfaceDimension);
	std::size_t spatialGridSize = std::accumulate(spatialGridLinearSize.begin(),
                                                  spatialGridLinearSize.begin() + this->surfaceDimension, 1.,
                                                  std::multiplies<>{});

	for (size_t i=0; i<this->allowedOrientations.size(); i++)
		delete this->voxels[i];
	this->allocateVoxels(this->allowedOrientations.size()*spatialGridSize);
	this->activeTopLevelVoxels = new bool[this->allowedOrientations.size()*fullSpatialGridSize];
	for(size_t i = 0; i < this->allowedOrientations.size()*fullSpatialGridSize; i++){
		this->activeTopLevelVoxels[i] = true;
	}
	this->spatialVoxelSize = this->initialVoxelSize;
	this->voxelNeighbourGrid = new NeighbourGrid<Voxel>(this->surfaceDimension,
                                                        this->spatialVoxelSize * fullSpatialGridLinearSize,
                                                        this->spatialVoxelSize);

	std::cout << " alocating " << (this->allowedOrientations.size()*spatialGridSize) << " voxels, now checking " << std::flush;

	size_t dotEvery = (spatialGridSize / 100) + 1;
	_OMP_PARALLEL_FOR
	for (size_t spatialIndex = 0; spatialIndex < spatialGridSize; spatialIndex++){
		RSAVector position;

		std::array<int, RSA_SPATIAL_DIMENSION> spatialGridVector{};

		calculateGrid<RSA_SPATIAL_DIMENSION>(spatialGridVector, spatialIndex, spatialGridLinearSize);
		for(unsigned char i=0; i<this->surfaceDimension; i++){
			position[i] = this->spatialVoxelSize * spatialGridVector[i]; // position point to the "left bottom" corner of a voxel
		}
		bool removeTopLevelVoxel = true;
		for(size_t angularIndex=0; angularIndex<this->allowedOrientations.size(); angularIndex++){
			size_t index = spatialIndex * this->allowedOrientations.size() + angularIndex;
			this->voxels[index] = new Voxel(position, this->allowedOrientations[angularIndex]);
			RSAOrientation aSize;
			for (unsigned short i=0; i<RSA_ANGULAR_DIMENSION; i++) {
				aSize[i]=0;
			}
			if (this->analyzeVoxel(this->voxels[index], nl, bc, this->spatialVoxelSize, aSize)){ // dividing only not overlapping voxels
				delete this->voxels[index];
				this->voxels[index] = nullptr;
			}else{
				removeTopLevelVoxel = false;
				this->voxels[index]->index = index;
			}
		}
		if (removeTopLevelVoxel == true)
			this->removeTopLevelVoxel(this->voxels[spatialIndex * this->allowedOrientations.size()]);
		if (spatialIndex%dotEvery == 0){ std::cout << "."; std::cout.flush(); }
	}

	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);
	//this->beginningVoxelNumber = ss*sa;
	this->length = this->allowedOrientations.size() * spatialGridSize;// this->beginningVoxelNumber;
	this->compactVoxelArray();

	this->voxelsInitialized = true;

	this->checkTopLevelVoxels();
	this->rebuildNeighbourGrid();
	return spatialGridSize * this->allowedOrientations.size();
}

Voxel *DiscreteAngleVoxelList::getVoxel(const RSAVector &pos, const RSAOrientation &angle) const{
	if (!this->voxelsInitialized){
		for(size_t i=0; i<this->length; i++){
			if(this->voxels[i]->getOrientation() == angle)
				return this->voxels[i];
		}
		return nullptr;
	}
	std::vector<Voxel *> *vTmp = this->voxelNeighbourGrid->getCell(pos);
	for(Voxel *v : *vTmp){
		if (v->isInside(pos, this->spatialVoxelSize, angle, this->angularVoxelSize)){
			return v;
		}
	}
	return nullptr;
}


void DiscreteAngleVoxelList::splitVoxel(Voxel *v, double spatialSize, const RSAOrientation &angularSize, Voxel **vRes) const
{
	unsigned short spatialLoop = 1 << this->surfaceDimension;

    RSAVector position = v->getPosition();
    RSAVector vpos = v->getPosition();
    RSAOrientation orientation = v->getOrientation();

    std::array<int, RSA_SPATIAL_DIMENSION> inpos{};
    inpos.fill(0);

	for(unsigned short i=0; i<spatialLoop; i++){
		for(unsigned short j=0; j < this->surfaceDimension; j++){
			position[j] = vpos[j] + inpos[j]*spatialSize;
		}
		vRes[i] = new Voxel(position, orientation);
		increment(inpos.data(), this->surfaceDimension, (unsigned char)1);
	} // for i
}

void DiscreteAngleVoxelList::getRandomPositionAndOrientation(RSAVector *position, RSAOrientation *orientation, Voxel *v, RND *rnd) const
{
	RSAVector vpos = v->getPosition();
	RSAOrientation vangle = v->getOrientation();

	for (unsigned short i=0; i < this->surfaceDimension; i++)
        (*position)[i] = vpos[i] + rnd->nextValue(this->spatialDistribution);
	for (unsigned short i=this->surfaceDimension; i < RSA_SPATIAL_DIMENSION; i++)
        (*position)[i] = 0.0;

	for (unsigned short i=0; i < RSA_ANGULAR_DIMENSION; i++)
        (*orientation)[i] = vangle[i];
}


RSAOrientation DiscreteAngleVoxelList::getAngularVoxelSize() const{
	RSAOrientation aSize;
	for (unsigned short i=0; i<RSA_ANGULAR_DIMENSION; i++)
		aSize[i] = 0;
	return aSize;
}

double DiscreteAngleVoxelList::getVoxelsVolume() const{
	double result = 0;
	for(size_t i = 0; i< this->length; i++){
		double s = 1.0;
		Voxel *v = this->voxels[i];
		RSAVector position = v->getPosition();
		for(unsigned short j=0; j<this->surfaceDimension; j++){
			if (position[j]+this->spatialVoxelSize > this->spatialRange){
				s *= this->spatialRange - position[j];
			}else{
				s *= this->spatialVoxelSize;
			}
		}
		result += s/static_cast<double>(this->allowedOrientations.size());
	}
	return result;
}
