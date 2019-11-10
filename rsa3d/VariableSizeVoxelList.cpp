/*
 * VariableSizeVoxelList.cpp
 *
 *  Created on: Nov 4, 2019
 *      Author: ciesla
 */

#include "VariableSizeVoxelList.h"
#include "utils/OMPMacros.h"

#include <math.h>


VariableSizeVoxelList::VariableSizeVoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange, double requestedAngularVoxelSize) : VoxelList(dim, packingSpatialSize, requestedSpatialVoxelSize, shapeAngularRange, requestedAngularVoxelSize){
	this->voxelMap = new double[1];
	this->voxelsDivisionCounters = new unsigned short[1];
	this->u01Distribution = new std::uniform_real_distribution<double>(0.0, 1.0);
}

VariableSizeVoxelList::~VariableSizeVoxelList() {
	delete this->voxelMap;
	delete this->voxelsDivisionCounters;
}

void VariableSizeVoxelList::allocateVoxels(size_t size){
	delete[] this->voxels;
	delete this->voxelsDivisionCounters;

	this->voxels = new Voxel*[size];
	this->voxelsDivisionCounters = new unsigned short[size];
	for(size_t i=0; i<size; i++){
		this->voxels[i] = nullptr;
		this->voxelsDivisionCounters[i] = 0;
	}
}

double VariableSizeVoxelList::getSpatialVoxelSize(size_t i){
	if (!this->voxelsInitialized)
		return this->spatialRange;
	else
		return this->initialVoxelSize/std::pow(2.0, this->voxelsDivisionCounters[i]);
}

double VariableSizeVoxelList::getAngularVoxelSize(size_t i){
	if (!this->voxelsInitialized)
		return this->angularRange;
	else
		return this->initialAngularVoxelSize/std::pow(2.0, this->voxelsDivisionCounters[i]);
}

double VariableSizeVoxelList::getVoxelVolume(size_t i){
	double spatialVolume = std::pow(this->getSpatialVoxelSize(i), this->surfaceDimension);
	double angularVolume = std::pow(this->getAngularVoxelSize(i), RSA_ANGULAR_DIMENSION);
	return (spatialVolume*angularVolume);
}

void VariableSizeVoxelList::createVoxelMap(){
	double sum = 0.0;
	delete this->voxelMap;
	this->voxelMap = new double[this->length];
	for(size_t i=0; i<this->length; i++){
		sum += this->getVoxelVolume(i);
		this->voxelMap[i] = sum;
	}
	for(size_t i=0; i<this->length; i++){
		this->voxelMap[i] /= sum;
	}
}

Voxel * VariableSizeVoxelList::getRandomVoxel(RND *rnd){
	double d = rnd->nextValue();
	size_t i1=0, i2=this->length-1, i;
	while(i2>i1+1 && i2>0){
		i = 0.5*(i1 + i2);
		if (this->voxelMap[i] > d)
			i2 = i-1;
		else
			i1 = i;
	}
	return this->voxels[i1];
}

void VariableSizeVoxelList::getRandomPositionAndOrientation(RSAVector *position, RSAOrientation *orientation, Voxel *v, RND *rnd){
	RSAVector vpos = v->getPosition();
	RSAOrientation vangle = v->getOrientation();

	for (unsigned short i=0; i < this->surfaceDimension; i++)
        (*position)[i] = vpos[i] + rnd->nextValue(this->u01Distribution)*this->getSpatialVoxelSize(v->index);
	for (unsigned short i=this->surfaceDimension; i < RSA_SPATIAL_DIMENSION; i++)
        (*position)[i] = 0.0;

	for (unsigned short i=0; i < RSA_ANGULAR_DIMENSION; i++)
        (*orientation)[i] = vangle[i] + rnd->nextValue(this->u01Distribution)*this->getAngularVoxelSize(v->index);
}


double VariableSizeVoxelList::splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc){
	if (this->disabled)
		return -1;

	size_t v0 = this->length;
	if (!this->voxelsInitialized){
		v0 = this->initVoxels(bc, nl);
		this->createVoxelMap();
		return ((double)this->length / (double)v0);
	}
	size_t voxelsFactor = (size_t)round( pow(2, this->surfaceDimension+RSA_ANGULAR_DIMENSION) );

	size_t newListSize = voxelsFactor*(this->length);
	Voxel** newList = new Voxel*[newListSize];
	unsigned short *newVoxelsDivisionCounters = new unsigned short[newListSize];
	for(size_t i=0; i<newListSize; i++){
		newList[i] = nullptr;
		newVoxelsDivisionCounters[i] = 0;
	}


	// temporary metrices for voxels after divisions. One separate matrix for each thread
	Voxel ***aVoxels = new Voxel**[_OMP_MAXTHREADS];
	size_t newVoxelsCounter[_OMP_MAXTHREADS];
	size_t initialNewVoxelsCounter[_OMP_MAXTHREADS];
	for(unsigned short i=0; i<_OMP_MAXTHREADS; i++){
		aVoxels[i] = new Voxel*[ voxelsFactor ];
		newVoxelsCounter[i] = 0;
		initialNewVoxelsCounter[i] = 0;
	}

	_OMP_PARALLEL_FOR
	for(size_t i=0; i<this->length; i++){
		// voxel is tested if it should remain active and if so it is divided
		if (!this->analyzeVoxel(this->voxels[i], nl, bc, this->getSpatialVoxelSize(i), this->getAngularVoxelSize(i))){
			// if too much new voxels there is no point in further splitting
			if (newVoxelsCounter[_OMP_THREAD_ID]*_OMP_MAXTHREADS <= maxVoxels){
				// preparing array of new voxels after division of this->voxels[i] (aVoxels)
				this->splitVoxel(this->voxels[i], this->getSpatialVoxelSize(i)/2.0, this->getAngularVoxelSize(i)/2.0, aVoxels[_OMP_THREAD_ID]);
				// analyzing new voxels, only the non covered ones will be added to the new array (newList)
				initialNewVoxelsCounter[_OMP_THREAD_ID] = newVoxelsCounter[_OMP_THREAD_ID];
				for(size_t j=0; j<voxelsFactor; j++){
					Voxel *v = aVoxels[_OMP_THREAD_ID][j];
					if( nl==nullptr || bc==nullptr || !this->analyzeVoxel(v, nl, bc, this->getSpatialVoxelSize(i)/2.0, this->getAngularVoxelSize(i)/2.0) ){
						if(this->voxels[i]->depth > 0){
							v->depth = this->voxels[i]->depth-1;
						}
						newList[i*voxelsFactor + j] = v;
						newVoxelsDivisionCounters[i*voxelsFactor + j] = this->voxelsDivisionCounters[i] + 1;
						newVoxelsCounter[_OMP_THREAD_ID]++;
					}else{
						// covered voxels are removed
						delete aVoxels[_OMP_THREAD_ID][j];
					}
				}
				if( (newVoxelsCounter[_OMP_THREAD_ID] - initialNewVoxelsCounter[_OMP_THREAD_ID]) > 0.25*voxelsFactor && newListSize > maxVoxels){
					// too many new voxels - the old one will remain
					for (size_t j=0; j<voxelsFactor; j++){
						if (newList[i*voxelsFactor + j]!=nullptr){
							delete newList[i*voxelsFactor + j];
							newList[i*voxelsFactor + j] = nullptr;
							newVoxelsDivisionCounters[i*voxelsFactor + j] = 0;
							newVoxelsCounter[_OMP_THREAD_ID]--;
						}
					}
					Assert(newVoxelsCounter[_OMP_THREAD_ID] == initialNewVoxelsCounter[_OMP_THREAD_ID]);
					newList[i*voxelsFactor + 0] = this->voxels[i];
					newVoxelsDivisionCounters[i*voxelsFactor + 0] = this->voxelsDivisionCounters[i];
					newVoxelsCounter[_OMP_THREAD_ID]++;
				}
			}else{ // no further splitting
				newList[i*voxelsFactor + 0] = this->voxels[i];
				newVoxelsDivisionCounters[i*voxelsFactor + 0] = this->voxelsDivisionCounters[i];
				newVoxelsCounter[_OMP_THREAD_ID]++;
			}
		}else{
			// original covered voxels are cleaned
			delete this->voxels[i];
			this->voxels[i] = nullptr;
		}
		if (i%10000 == 0){ std::cout << "." << std::flush; }
	}

	// delete temporary thread matrices. Covered voxels have been already removed
	for(int i=0; i<_OMP_MAXTHREADS; i++){
		delete[] aVoxels[i];
	}
	delete[] aVoxels;

	// clearing the old list and processing the new one
	for(size_t i=0; i<this->length; i++){
		if(this->voxels[i]!=nullptr && this->voxels[i]!=newList[i*voxelsFactor + 0])
			delete this->voxels[i];
	}
	delete[] this->voxels;
	this->voxels = newList;
	delete this->voxelsDivisionCounters;
	this->voxelsDivisionCounters = newVoxelsDivisionCounters;
	this->length = newListSize;
	size_t newVoxelsSize = 0;
	for(unsigned short i=0; i<_OMP_MAXTHREADS; i++)
		newVoxelsSize += newVoxelsCounter[i];

	std::cout << " compacting" << std::flush;
	this->length = this->compactVoxelArray();
	Assert(this->length == newVoxelsSize);

	this->createVoxelMap();
	this->rebuildNeighbourGrid();
	//	this->checkTopLevelVoxels();
	return ((double)this->length / (double)v0);
}

void VariableSizeVoxelList::moveVoxelInList(size_t from, size_t to){
	Assert(this->voxels[to] == nullptr);
	this->voxels[to] = this->voxels[from];
	this->voxels[from] = nullptr;
	this->voxelsDivisionCounters[to] = this->voxelsDivisionCounters[from];
	this->voxelsDivisionCounters[from] = 0;
}

double VariableSizeVoxelList::getVoxelsVolume(){
	double result = 0;
	for(size_t i = 0; i< this->length; i++){
		double s = 1.0;
		Voxel *v = this->voxels[i];
		RSAVector position = v->getPosition();
		for(unsigned short j=0; j<this->surfaceDimension; j++){
			if (position[j]+this->getSpatialVoxelSize(i) > this->spatialRange){
				s *= this->spatialRange - position[j];
			}else{
				s *= this->getSpatialVoxelSize(i);
			}
		}
	    double a = 1.0;
		RSAOrientation orientation = v->getOrientation();
		for(unsigned char j=0; j<RSA_ANGULAR_DIMENSION; j++){
			if (orientation[j]+this->getAngularVoxelSize(i) > this->angularRange){
				a *= this->angularRange - orientation[j];
			}else{
				a *= this->getAngularVoxelSize(i);
			}
			a /= this->angularRange;
		}
		result += s*a;
	}
	return result;
}

Voxel *VariableSizeVoxelList::getVoxel(const RSAVector &pos, const RSAOrientation &angle){
	if (!this->voxelsInitialized)
		return this->voxels[0];
	std::vector<Voxel *> *vTmp = this->voxelNeighbourGrid->getCell(pos);
	for(Voxel *v : *vTmp){
		if (v->isInside(pos, this->getSpatialVoxelSize(v->index), angle, this->getAngularVoxelSize(v->index))){
			return v;
		}
	}
	return nullptr;
}

