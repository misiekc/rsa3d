/*
 * VoxelList.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#include <iostream>
#include <cmath>
#include <algorithm>

#include "VoxelList.h"
#include "utils/OMPMacros.h"
#include "utils/Assertions.h"

/**
 * d - requested initial size of a voxel
 */

double VoxelList::findFloorSize(double d){
	double dRes = 1.0;
	while (dRes > d)
		dRes /= 2;
	while (2*dRes < d)
		dRes *= 2;
	return dRes;
}


double VoxelList::findCeilSize(double d){
	double dRes = 1.0;
	while (dRes < d)
		dRes *= 2;
	while (dRes > d)
		dRes /= 2;
	return 2*dRes;
}

size_t VoxelList::findArraySize(double range, double cellSize) const {
	return (size_t)(range/cellSize) + 1;
}

void VoxelList::allocateVoxels(size_t size){
	delete[] this->voxels;
	this->voxels = new Voxel*[size];
	for(size_t i=0; i<size; i++){
		this->voxels[i] = nullptr;
	}
}

/**
 * Dummy constructor for derived classes
 */
VoxelList::VoxelList(){
	this->surfaceDimension = RSA_SPATIAL_DIMENSION;
	this->spatialRange = 0.0;
	this->angularRange = 0.0;
	this->angularVoxelSize = this->angularRange;


	this->initialVoxelSize = 0.0;
	this->initialAngularVoxelSize = 0.0;
	this->spatialVoxelSize = 0.0;
	//this->beginningVoxelNumber = 1;
	this->activeTopLevelVoxels = nullptr;
	this->voxelNeighbourGrid = nullptr;


	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);
	this->disabled = true;
	this->voxelsInitialized = false;

	this->voxels = nullptr;
	this->length = 0;
}



VoxelList::VoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange, double requestedAngularVoxelSize){
	Expects(dim > 0);
	Expects(dim <= RSA_SPATIAL_DIMENSION);
	Expects(packingSpatialSize > 0.0);
	Expects(requestedSpatialVoxelSize > 0.0);
	Expects(shapeAngularRange > 0.0);
	Expects(requestedAngularVoxelSize > 0.0);

	this->surfaceDimension = dim;
	this->spatialRange = packingSpatialSize;
	this->angularRange = shapeAngularRange;
	this->angularVoxelSize = this->angularRange;


	this->initialVoxelSize = this->findFloorSize(requestedSpatialVoxelSize);
	this->initialAngularVoxelSize = this->findCeilSize(requestedAngularVoxelSize);
	this->spatialVoxelSize = this->spatialRange;
	//this->beginningVoxelNumber = 1;
	this->activeTopLevelVoxels = nullptr;
	this->voxelNeighbourGrid = nullptr;


	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);
	this->disabled = false;
	this->voxelsInitialized = false;

	RSAVector position;
	RSAOrientation orientation;
	for(unsigned short i=0; i<RSA_SPATIAL_DIMENSION; i++)
		position[i] = 0.0;
	for(unsigned short i=0; i<RSA_ANGULAR_DIMENSION; i++)
		orientation[i] = 0.0;

	this->voxels = new Voxel*[1];
	this->voxels[0] = new Voxel(position, orientation);
	this->voxels[0]->index = 0;
	this->length = 1;
}

VoxelList::~VoxelList() {
	delete this->spatialDistribution;
	delete this->angularDistribution;

	for(size_t i=0; i<this->length; i++){
		delete this->voxels[i];
	}
	delete[] this->voxels;
	delete[] this->activeTopLevelVoxels;
	delete this->voxelNeighbourGrid;
}

void VoxelList::disable(){
	this->disabled = true;
}

unsigned int VoxelList::initVoxels(RSABoundaryConditions *bc, NeighbourGrid<const RSAShape> *nl){
    // "Full spatial" variables correspond to the whole cube simulation box (for initial voxels and neighbour grid)
    // "Spatial" (without "full") variables correspond to the acutal part where voxels will be created (for ordinary
    // voxels)
	size_t fullSpatialGridLinearSize = this->findArraySize(this->spatialRange, this->initialVoxelSize);
	size_t angularGridLinearSize = this->findArraySize(this->angularRange, this->initialAngularVoxelSize);
	size_t fullSpatialGridSize = (size_t)(pow(fullSpatialGridLinearSize, this->surfaceDimension) + 0.5);
	size_t angularGridSize = (size_t)(pow(angularGridLinearSize, RSA_ANGULAR_DIMENSION) + 0.5);
	std::array<std::size_t, RSA_SPATIAL_DIMENSION> spatialGridLinearSize = this->calculateSpatialGridLinearSize();
	Assert(RSA_SPATIAL_DIMENSION >= this->surfaceDimension);
	std::size_t spatialGridSize = std::accumulate(spatialGridLinearSize.begin(),
                                                  spatialGridLinearSize.begin() + this->surfaceDimension, 1.,
                                                  std::multiplies<>{});

	delete this->voxels[0];
	this->allocateVoxels(spatialGridSize * angularGridSize);
	this->activeTopLevelVoxels = new bool[fullSpatialGridSize];
	for(size_t i = 0; i < fullSpatialGridSize; i++){
		this->activeTopLevelVoxels[i] = true;
	}
	this->spatialVoxelSize = this->initialVoxelSize;
	this->angularVoxelSize = this->initialAngularVoxelSize;
	this->voxelNeighbourGrid = new NeighbourGrid<Voxel>(this->surfaceDimension,
                                                        this->spatialVoxelSize * fullSpatialGridLinearSize,
                                                        this->spatialVoxelSize);

	std::cout << " alocating " << (spatialGridSize * angularGridSize) << " voxels, now checking " << std::flush;

	size_t dotEvery = (spatialGridSize / 100) + 1;
	_OMP_PARALLEL_FOR
	for (size_t spatialIndex = 0; spatialIndex < spatialGridSize; spatialIndex++){
		RSAVector position;
		RSAOrientation orientation{};

		std::array<int, RSA_SPATIAL_DIMENSION> spatialGridVector{};
		std::array<int, RSA_ANGULAR_DIMENSION> angularGridVector{};

		calculateGrid<RSA_SPATIAL_DIMENSION>(spatialGridVector, spatialIndex, spatialGridLinearSize);
		for(unsigned char i=0; i<this->surfaceDimension; i++){
			position[i] = this->spatialVoxelSize * spatialGridVector[i]; // position point to the "left bottom" corner of a voxel
		}
		bool removeTopLevelVoxel = true;
		for(size_t angularIndex = 0; angularIndex < angularGridSize; angularIndex++ ){
			calculateGrid<RSA_ANGULAR_DIMENSION>(angularGridVector, angularIndex, angularGridLinearSize);
			for(unsigned char i=0; i<RSA_ANGULAR_DIMENSION; i++){
				orientation[i] = this->angularVoxelSize * angularGridVector[i]; // orientation point to the "left bottom" corner of a voxel
			}
			size_t index = spatialIndex * angularGridSize + angularIndex;
			this->voxels[index] = new Voxel(position, orientation);
			if (this->analyzeVoxel(this->voxels[index], nl, bc, this->spatialVoxelSize, this->angularVoxelSize)){ // dividing only not overlapping voxels
				delete this->voxels[index];
				this->voxels[index] = nullptr;
			}else{
				removeTopLevelVoxel = false;
				this->voxels[index]->index = index;
			}
		}
		if (removeTopLevelVoxel == true)
			this->removeTopLevelVoxel(this->voxels[spatialIndex * angularGridSize]);
		if (spatialIndex%dotEvery == 0){ std::cout << "."; std::cout.flush(); }
	}

	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);
	delete this->angularDistribution;
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);

	//this->beginningVoxelNumber = ss*sa;
	this->length = angularGridSize * spatialGridSize;// this->beginningVoxelNumber;
	this->compactVoxelArray();

	this->voxelsInitialized = true;

	this->checkTopLevelVoxels();
	this->rebuildNeighbourGrid();
	return spatialGridSize * angularGridSize;
}

void VoxelList::rebuildNeighbourGrid(){
	this->voxelNeighbourGrid->clear();
	for(size_t i=0; i<this->length; i++){
		this->voxelNeighbourGrid->add(this->voxels[i], this->voxels[i]->getPosition());
	}
}


void VoxelList::getNeighbours(std::vector<Voxel *> *result, const RSAVector &da){
	return this->voxelNeighbourGrid->getNeighbours(result, da);
}

void VoxelList::moveVoxelInList(size_t from, size_t to){
	Assert(this->voxels[to] == nullptr);
	this->voxels[to] = this->voxels[from];
	this->voxels[from] = nullptr;
}

void VoxelList::compactVoxelArray(){
	if (this->length==0)
		return;

	size_t counter = 0;
	for(size_t i=0; i<this->length; i++)
		if (this->voxels[i]!=nullptr)
			counter++;

	if (counter==0 || counter==this->length){
		this->length = counter;
		return;
	}

	size_t beginIndex = 0;
	size_t endIndex = this->length-1;

	do{
		// moving lastIndex to the last non-null value
		while(endIndex > 0 && this->voxels[endIndex] == nullptr)
			(endIndex)--;

		// moving beginIdex to the first null value
		while(beginIndex < this->length && this->voxels[beginIndex] !=nullptr)
			beginIndex++;

		if(beginIndex<endIndex){
			this->moveVoxelInList(endIndex, beginIndex);
			beginIndex++;
			endIndex--;
		}
	}while(beginIndex<=counter && endIndex>=counter);

	Assert(counter == endIndex+1);

	for(size_t i=0; i<=endIndex; i++)
		this->voxels[i]->index = i;

	this->length = endIndex+1;
	this->refreshTopLevelVoxels();
}

int VoxelList::getIndexOfTopLevelVoxel(const RSAVector &da) const{

	double daArray[RSA_SPATIAL_DIMENSION];
	da.copyToArray(daArray);
	int n = (int)(this->spatialRange/this->initialVoxelSize) + 1;
	int index = position2i(daArray, this->surfaceDimension, n*this->initialVoxelSize, this->initialVoxelSize, n);
	return index;
}


void VoxelList::removeTopLevelVoxel(Voxel *v){
	if (this->disabled || !this->voxelsInitialized)
		return;
	this->activeTopLevelVoxels[this->getIndexOfTopLevelVoxel(v->getPosition())]=false;
}


void VoxelList::checkTopLevelVoxels(){

	RND rnd;
	RSAVector pos;
	RSAOrientation angle{};

	for(size_t i=0; i<this->length; i++){
		Voxel *v = this->voxels[i];
		this->getRandomPositionAndOrientation(&pos, &angle, v, &rnd);
		int i1 = this->getIndexOfTopLevelVoxel(pos);
		int i2 = this->getIndexOfTopLevelVoxel(v->getPosition());
		if (i1!=i2){
			std::cout << "checkTopVoxels problem" << std::endl;
		}
		for(int j=0; j<10; j++){
			this->getRandomPositionAndOrientation(&pos, &angle, v, &rnd);
			i2 = this->getIndexOfTopLevelVoxel(pos);
			if (i1!=i2){
				std::cout << "checkTopVoxels problem" << std::endl;
			}
		}

	}
}

size_t VoxelList::countActiveTopLevelVoxels(){
	size_t result = 0;
	size_t n = (size_t)(this->spatialRange/this->initialVoxelSize) + 1;
	auto max = static_cast<std::size_t>(std::round(std::pow(n, this->surfaceDimension)));
	for(size_t i=0; i<max; i++){
		if(this->activeTopLevelVoxels[i])
			result++;
	}
	return result;
}

void VoxelList::refreshTopLevelVoxels(){
	size_t n = (size_t)(this->spatialRange/this->initialVoxelSize) + 1;
    auto max = static_cast<std::size_t>(std::round(std::pow(n, this->surfaceDimension)));
	for(size_t i=0; i<max; i++){
		this->activeTopLevelVoxels[i] = false;
	}

	_OMP_PARALLEL_FOR
	for(size_t i=0; i<this->length; i++){
		size_t index = this->getIndexOfTopLevelVoxel(this->voxels[i]->getPosition());
		this->activeTopLevelVoxels[index] = true;
	}

}


Voxel *VoxelList::getVoxel(const RSAVector &pos, const RSAOrientation &angle) const{
	if (!this->voxelsInitialized)
		return this->voxels[0];
	std::vector<Voxel *> *vTmp = this->voxelNeighbourGrid->getCell(pos);
	for(Voxel *v : *vTmp){
		if (v->isInside(pos, this->spatialVoxelSize, angle, this->angularVoxelSize)){
			return v;
		}
	}
	return nullptr;
}

void VoxelList::splitVoxel(Voxel *v, double spatialSize, double angularSize, Voxel **vRes) const{

	unsigned short spatialLoop = 1 << this->surfaceDimension;
	unsigned short angularLoop = 1 << RSA_ANGULAR_DIMENSION;

    RSAVector position = v->getPosition();
    RSAVector vpos = v->getPosition();
    RSAOrientation orientation = v->getOrientation();
    RSAOrientation vangle = v->getOrientation();

    std::array<int, RSA_SPATIAL_DIMENSION> inpos{};
    std::array<int, RSA_ANGULAR_DIMENSION> inangle{};
    inpos.fill(0);
    inangle.fill(0);

	for(unsigned short i=0; i<spatialLoop; i++){
		for(unsigned short j=0; j < this->surfaceDimension; j++){
			position[j] = vpos[j] + inpos[j]*spatialSize;
		}
		inangle.fill(0);
		for(unsigned short j=0; j<angularLoop; j++){
			for(unsigned short k=0; k<RSA_ANGULAR_DIMENSION; k++){
				orientation[k] = vangle[k] + inangle[k]*angularSize;
			}
			vRes[i*angularLoop + j] = new Voxel(position, orientation);
			increment(inangle.data(), RSA_ANGULAR_DIMENSION, (unsigned char)1);
		} // for j
		increment(inpos.data(), this->surfaceDimension, (unsigned char)1);
	} // for i
}


bool VoxelList::isVoxelInsidePacking(const Voxel *v, [[maybe_unused]] double spatialSize) const {
	RSAVector vpos = v->getPosition();
	RSAOrientation vangle = v->getOrientation();
	for(unsigned short i=0; i < this->surfaceDimension; i++){
		if (vpos[i] >= this->spatialRange){
			return false;
		}
	}
	for(unsigned short i=0; i<RSA_ANGULAR_DIMENSION; i++){
		if (vangle[i] >= this->angularRange){
			return false;
		}
	}
	return true;
}

/**
 * returns true when the whole voxel is inside an exclusion area of any shape in shapes
 * To determine it the method tires to split voxel up to level of maxDepth
 */
bool VoxelList::isVoxelInsideExclusionZone(Voxel *v, double spatialSize, double angularSize,
										   std::vector<const RSAShape *> *shapes, RSABoundaryConditions *bc,
                                           unsigned short depth) const{

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
			if (finalVoxels[j]!=v)
				delete finalVoxels[j]; //deleting virtual voxel which were already splitted
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
		for(const RSAShape *s : *shapes){
			// jesli choc jeden z mniejszych jest false zwracamy false
			isInside = s->voxelInside(bc, finalVoxels[i]->getPosition(), finalVoxels[i]->getOrientation(), ss, as);
			if (isInside)
				break; // znaleziony, przechodzimy do nastepnego woksela
		}
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

// old recursive version
// leaved for testing purposes
bool VoxelList::isVoxelInsideExclusionZoneOld(Voxel *v, double spatialSize, double angularSize,
										   std::vector<const RSAShape *> *shapes, RSABoundaryConditions *bc,
                                           unsigned short depth){
	// if voxel is outside the packing it is inside exclusion zone
//	if (!this->isVoxelInsidePacking(v))
//		return true;
	// otherwise checking

	bool isInside = false;
	for(const RSAShape *s : *shapes){
		isInside = s->voxelInside(bc, v->getPosition(), v->getOrientation(), spatialSize, angularSize);
		if (isInside)
			break;
	}

	if (isInside || depth == 0){
		return isInside;
	// if cannot determine that it is inside and depth > 0 split and recursively check children
	}else{
		// dzielimy voxel na mniejsze
		double ss = spatialSize / 2.0;
		double as = angularSize / 2.0;
		int arrayLenght = (int)round( pow(2, this->surfaceDimension+RSA_ANGULAR_DIMENSION) );
		Voxel **aVoxels = new Voxel*[ arrayLenght ];
		this->splitVoxel(v, ss, as, aVoxels);
		bool bRes = true;
		// sprawdzamy kazdy z mniejszych
		for(int i=0; i<arrayLenght; i++){
			// jesli choc jeden z mniejszych jest false zwracamy false
			if (!this->isVoxelInsideExclusionZone(aVoxels[i], ss, as, shapes, bc, depth-1)){
				bRes = false;
				break;
			}
		}
		// w przeciwnym razie zwracamy true;

		for(int i=0; i<arrayLenght; i++){
			delete aVoxels[i];
		}
		delete[] aVoxels;

		return bRes;
	}
}

bool VoxelList::isTopLevelVoxelActive(Voxel *v) const{

	// checking if initial voxel containing v is active (do not have a shape inside)
	if (!this->voxelsInitialized)
		return true;
	int index = this->getIndexOfTopLevelVoxel(v->getPosition());
	return this->activeTopLevelVoxels[index];
}

bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc, double spatialSize, double angularSize, unsigned short depth) const{
	if (!this->disabled){ // && (depth > v->depth || depth==0) ){

	    if (!isTopLevelVoxelActive(v) || !this->isVoxelInsidePacking(v, spatialSize) )
			return true;

	    std::vector<const RSAShape*> tmpShapes, shapes;
		    nl->getNeighbours(&tmpShapes, v->getPosition());

		int maxNo = v->lastAnalyzed;
		for(const RSAShape *s: tmpShapes){
			if (v->lastAnalyzed < s->no || depth > v->depth){
				shapes.push_back(s);
				if(maxNo < s->no)
					maxNo = s->no;
			}
		}
//		bool isInside = this->isVoxelInsideExclusionZoneOld(v, spatialSize, angularSize, &shapes, bc, depth);
		bool isInside = this->isVoxelInsideExclusionZone(v, spatialSize, angularSize, &shapes, bc, depth);
		v->depth = depth;
		v->lastAnalyzed = maxNo;
		return isInside;
	}
	return false;
}

size_t VoxelList::analyzeVoxels(RSABoundaryConditions *bc, NeighbourGrid<const RSAShape> *nl, unsigned short depth) {

	size_t begin = this->length;
	size_t dotEvery = (this->length/100)+1;


	_OMP_PARALLEL_FOR
	for (size_t i = 0; i < this->length; i++) {
		Voxel *v = this->voxels[i];
		bool bRemove = this->analyzeVoxel(v, nl, bc, this->spatialVoxelSize, this->angularVoxelSize, depth);
		if (bRemove){
			this->voxelNeighbourGrid->remove(v, v->getPosition());
			delete v;
			this->voxels[i] = nullptr;
		}
		if (i%dotEvery == 0){ std::cout << "."; std::cout.flush(); }
	}

	std::cout << " compacting" << std::flush;

	this->compactVoxelArray();

	return begin - this->length;
}

unsigned short VoxelList::splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc){
	if (this->disabled || this->spatialVoxelSize<2*minDx)
		return VoxelList::NO_SPLIT;

	if (!this->voxelsInitialized){
		this->initVoxels(bc, nl);
		return VoxelList::NO_SPLIT_BUT_INITIALIZED;
	}


	// number of created voxels
	size_t newVoxelsCounter = 0;

	size_t voxelsFactor = 1;
	if (this->spatialVoxelSize>0)
		voxelsFactor *= (size_t)round( pow(2, this->surfaceDimension));
	if (this->angularVoxelSize>0)
		voxelsFactor *= (size_t)round( pow(2, RSA_ANGULAR_DIMENSION));

	size_t newListSize = voxelsFactor*(this->length);
	Voxel** newList = new Voxel*[ newListSize ];
	for(size_t i=0; i<newListSize; i++){
		newList[i] = nullptr;
	}

	// temporary metrices for voxels after divisions. One separate matrix for each thread
	Voxel ***aVoxels = new Voxel**[_OMP_MAXTHREADS];
	for(int i=0; i<_OMP_MAXTHREADS; i++){
		aVoxels[i] = new Voxel*[ voxelsFactor ];
	}

	size_t dotEvery = (this->length/100)+1;

	_OMP_PARALLEL_FOR
	for(size_t i=0; i<this->length; i++){
		// voxel is tested if it should remain active and if so it is divided
		if (!this->analyzeVoxel(this->voxels[i], nl, bc, this->spatialVoxelSize, this->angularVoxelSize)){
			// if too much new voxels there is no point in further splitting
			if (newVoxelsCounter <= maxVoxels){
				// preparing array of new voxels after division of this->voxels[i] (aVoxels)
				this->splitVoxel(this->voxels[i], this->spatialVoxelSize/2.0, this->angularVoxelSize/2.0, aVoxels[_OMP_THREAD_ID]);
				// analyzing new voxels, only the non covered ones will be added to the new array (newList)
				for(size_t j=0; j<voxelsFactor; j++){
					Voxel *v = aVoxels[_OMP_THREAD_ID][j];
					if( nl==nullptr || bc==nullptr || !this->analyzeVoxel(v, nl, bc, this->spatialVoxelSize/2.0, this->angularVoxelSize/2.0) ){
						if(this->voxels[i]->depth > 0){
							v->depth = this->voxels[i]->depth-1;
						}
						newList[i*voxelsFactor + j] = v;
						_OMP_ATOMIC
						newVoxelsCounter++;
					}else{
						// covered voxels are removed
						delete aVoxels[_OMP_THREAD_ID][j];
					}
				} //for
			} // if
		}else{
			// original covered voxels are cleaned
			delete this->voxels[i];
			this->voxels[i] = nullptr;
		}
		if (i%dotEvery == 0){ std::cout << "." << std::flush; }
	}

	// delete temporary thread matrices. Covered voxels have been already removed
	for(int i=0; i<_OMP_MAXTHREADS; i++){
		delete[] aVoxels[i];
	}
	delete[] aVoxels;

	if (newVoxelsCounter > maxVoxels){
		// too much voxels. Voxel splitting cancelled - using oryginal list instead of the new one
		std::cout << " too many new voxels (" << newVoxelsCounter << ">" << maxVoxels <<"): - cancel splitting," << std::flush;
		for(size_t i=0; i<newListSize; i++){
			if (newList[i]!=nullptr)
				delete newList[i];
		}
		delete[] newList;

		std::cout << " compacting" << std::flush;
		this->compactVoxelArray();
		this->rebuildNeighbourGrid();
		return VoxelList::NO_SPLIT_DUE_TO_VOXELS_LIMIT;

	}else{
		// clearing the old list and processing the new one
		for(size_t i=0; i<this->length; i++){
			if(this->voxels[i]!=nullptr)
				delete this->voxels[i];
		}
		delete[] this->voxels;
		this->voxels = newList;
		this->length = newListSize;

		this->spatialVoxelSize = (this->spatialVoxelSize/2.0)*this->dxFactor;
		delete this->spatialDistribution;
		this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);

		this->angularVoxelSize = (this->angularVoxelSize/2.0)*this->dxFactor;
		delete this->angularDistribution;
		this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);


		std::cout << " compacting" << std::flush;
		this->compactVoxelArray();
		this->rebuildNeighbourGrid();
		return VoxelList::NORMAL_SPLIT;

	}
	//	this->checkTopLevelVoxels();
}

void VoxelList::getRandomEntry(RSAVector *position, RSAOrientation *orientation, Voxel **v, RND *rnd) {
	double d = rnd->nextValue();
	*v = (this->voxels[(int)(d*(this->length))]);
	this->getRandomPositionAndOrientation(position, orientation, *v, rnd);
}

Voxel * VoxelList::getVoxel(int i) {
	return this->voxels[i];
}

void VoxelList::getRandomPositionAndOrientation(RSAVector *position, RSAOrientation *orientation, Voxel *v, RND *rnd)
{
	RSAVector vpos = v->getPosition();
	RSAOrientation vangle = v->getOrientation();

	for (unsigned short i=0; i < this->surfaceDimension; i++)
        (*position)[i] = vpos[i] + rnd->nextValue(this->spatialDistribution);
	for (unsigned short i=this->surfaceDimension; i < RSA_SPATIAL_DIMENSION; i++)
        (*position)[i] = 0.0;

	for (unsigned short i=0; i < RSA_ANGULAR_DIMENSION; i++)
        (*orientation)[i] = vangle[i] + rnd->nextValue(this->angularDistribution);
}


double VoxelList::getSpatialVoxelSize() const {
	return this->spatialVoxelSize;
}


double VoxelList::getAngularVoxelSize() const {
	return this->angularVoxelSize;
}


Voxel* VoxelList::get(int i){
	return this->voxels[i];
}


size_t VoxelList::getLength() const{
	return this->length;
}


double VoxelList::getVoxelsVolume() const{
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
	    double a = 1.0;
		RSAOrientation orientation = v->getOrientation();
		for(unsigned char j=0; j<RSA_ANGULAR_DIMENSION; j++){
			if (orientation[j]+this->angularVoxelSize > this->angularRange){
				a *= this->angularRange - orientation[j];
			}else{
				a *= this->angularVoxelSize;
			}
			a /= this->angularRange;
		}
		result += s*a;
	}
	return result;
}


std::string VoxelList::toPovray() const{
	std::string sRes = "";

	for(size_t i=0; i<this->length; i++){
		sRes += this->voxels[i]->toPovray(this->spatialVoxelSize) + "\r\n";
	}
	return sRes;
}


std::string VoxelList::toWolfram() const{
	std::stringstream out;

	for(size_t i=0; i<this->length; i++){
		out << this->voxels[i]->toWolfram(this->spatialVoxelSize, this->angularVoxelSize);
		if (i!=this->length-1)
			out << ", ";
		out << std::endl;
		if (RSA_ANGULAR_DIMENSION > 0){
			out << "(* angles: [ " << this->voxels[i]->getOrientation()[0] << ", " << (this->voxels[i]->getOrientation()[0] + this->angularVoxelSize) << ") *)" << std::endl;
		}
	}
	return out.str();
}


void VoxelList::store(std::ostream &f) const{
	unsigned short sd = RSA_SPATIAL_DIMENSION;
	unsigned short ad = RSA_ANGULAR_DIMENSION;
	f.write((char *)(&sd), sizeof(unsigned char));
	if (ad>0)
		f.write((char *)(&ad), sizeof(unsigned char));
	f.write((char *) &this->spatialVoxelSize, sizeof(double));
	f.write((char *) &this->angularVoxelSize, sizeof(double));

	size_t size = this->length;
	f.write((char *)(&size), sizeof(size_t));
	for(size_t i=0; i<size; i++){
		this->voxels[i]->store(f);
	}
}


void VoxelList::restore(std::istream &f){
	unsigned char sd = RSA_SPATIAL_DIMENSION;
	unsigned char ad = RSA_ANGULAR_DIMENSION;

	f.read((char *)(&sd), sizeof(unsigned char));
	if (ad > 0)
		f.read((char *)(&ad), sizeof(unsigned char));

	if (sd!=RSA_SPATIAL_DIMENSION || ad!=RSA_ANGULAR_DIMENSION){
		std::cout << "[ERROR] cannot restore VoxelList: incompatible dimensions: read " << f.gcount() << " bytes." << std::endl;
		return;
	}
	f.read((char *) &(this->spatialVoxelSize), sizeof(double));
	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);

	f.read((char *)&this->angularVoxelSize, sizeof(double));
	delete this->angularDistribution;
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);

	size_t ns = this->findArraySize(this->spatialRange, this->initialVoxelSize);
	this->voxelNeighbourGrid = new NeighbourGrid<Voxel>(this->surfaceDimension, this->initialVoxelSize*ns, this->initialVoxelSize);
	this->voxelNeighbourGrid->clear();
	size_t ss = (size_t)(pow(ns, this->surfaceDimension)+0.5);

	this->activeTopLevelVoxels = new bool[ss];
	for(size_t i = 0; i<ss; i++){
		this->activeTopLevelVoxels[i] = false;
	}

	for(size_t i=0; i<this->length; i++){
		delete this->voxels[i];
	}
	delete[] this->voxels;

	size_t size;
	int topIndex;
	f.read((char *)&size, sizeof(size_t));
	this->voxels = new Voxel *[size];
	for(size_t i=0; i<size; i++){
		Voxel *v = new Voxel();
		v->restore(f);
		this->voxels[i] = v;
		topIndex = this->getIndexOfTopLevelVoxel(v->getPosition());
		this->activeTopLevelVoxels[topIndex] = true;
	}
	this->length = size;
	this->compactVoxelArray();
	this->rebuildNeighbourGrid();
	this->voxelsInitialized = true;
}

std::array<std::size_t, RSA_SPATIAL_DIMENSION> VoxelList::calculateSpatialGridLinearSize() const {
    std::array<std::size_t, RSA_SPATIAL_DIMENSION> result{};
    result.fill(this->findArraySize(this->spatialRange, this->initialVoxelSize));
    return result;
}
