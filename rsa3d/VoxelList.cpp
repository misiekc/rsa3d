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

#ifdef _OPENMP
#include <omp.h>
#endif


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



inline size_t VoxelList::findArraySize(double range, double cellSize){
	return (size_t)(range/cellSize) + 1;
}

VoxelList::VoxelList(double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange, double requestedAngularVoxelSize){
	this->spatialRange = packingSpatialSize;
	this->spatialVoxelSize = this->findFloorSize(requestedSpatialVoxelSize);
	this->initialVoxelSize = this->spatialVoxelSize;

	this->angularRange = shapeAngularRange;
	this->angularVoxelSize = this->findCeilSize(requestedAngularVoxelSize);
	this->initialAngularVoxelSize = this->angularVoxelSize;

	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);
	this->disabled = false;

	size_t ns = this->findArraySize(this->spatialRange, this->spatialVoxelSize);
	size_t na = this->findArraySize(this->angularRange, this->angularVoxelSize);
	size_t ss = (size_t)(pow(ns, RSA_SPATIAL_DIMENSION)+0.5);
	size_t sa = (size_t)(pow(na, RSA_ANGULAR_DIMENSION)+0.5);
	this->voxels = new Voxel*[ss*sa];
	this->activeTopLevelVoxels = new bool[ss];
	this->voxelNeighbourGrid = new NeighbourGrid<Voxel>(this->spatialVoxelSize*ns, ns);

	this->initVoxels();
	for(size_t i = 0; i<ss; i++){
		this->activeTopLevelVoxels[i] = true;
	}


	this->spatialVoxelSize *= this->dxFactor;
	this->beginningVoxelNumber = ss*sa;
	this->length = this->beginningVoxelNumber;
	this->rebuildNeighbourGrid();
	}


VoxelList::~VoxelList() {
	delete this->spatialDistribution;

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

void VoxelList::initVoxels(){
	size_t ns = this->findArraySize(this->spatialRange, this->spatialVoxelSize);
	size_t na = this->findArraySize(this->angularRange, this->angularVoxelSize);
	RSAVector position;
	RSAOrientation orientation{};
	std::array<int, RSA_SPATIAL_DIMENSION> ins{};
	std::array<int, RSA_ANGULAR_DIMENSION> ina{};

    ins.fill(0);
	int index = 0;
	do{
		for(unsigned char i=0; i<RSA_SPATIAL_DIMENSION; i++){
			position[i] = this->spatialVoxelSize*ins[i]; // position point to the "left bottom" corner of a voxel
		}
/*
		int i = position2i(position, RSA_SPATIAL_DIMENSION, n*this->voxelSize, this->voxelSize, n);
		if(index!=i){
			std::cout << "VoxelList::initVoxels: Problem: " << index << " != " << i << std::endl;
		}
*/
		ina.fill(0);
		do{
			for(unsigned char i=0; i<RSA_ANGULAR_DIMENSION; i++){
				orientation[i] = this->angularVoxelSize*ina[i]; // orientation point to the "left bottom" corner of a voxel
			}
			this->voxels[index] = new Voxel(position, orientation);
			index++;
		}while(increment(ina.data(), RSA_ANGULAR_DIMENSION, na-1));
	}while(increment(ins.data(), RSA_SPATIAL_DIMENSION, ns-1));
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

void VoxelList::compactVoxelArray(Voxel **list, int &endIndex){
	int beginIndex = 0;

	while(list[endIndex] == nullptr && endIndex > -1)
		(endIndex)--;

	while (beginIndex<endIndex){
		if(list[beginIndex] == nullptr){
			list[beginIndex] = list[endIndex];
			list[endIndex] = nullptr;
		}
		beginIndex++;
		while(list[endIndex] == nullptr && endIndex > -1)
			endIndex--;
	}
}

int VoxelList::getIndexOfTopLevelVoxel(const RSAVector &da){

	double daArray[RSA_SPATIAL_DIMENSION];
	da.copyToArray(daArray);
	int n = (int)(this->spatialRange/this->initialVoxelSize) + 1;
	int index = position2i(daArray, RSA_SPATIAL_DIMENSION, n*this->initialVoxelSize, this->initialVoxelSize, n);
	return index;
}


void VoxelList::removeTopLevelVoxel(Voxel *v){
	if (this->disabled)
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
		int index = this->getIndexOfTopLevelVoxel(pos);
		if (this->getIndexOfTopLevelVoxel(v->getPosition())!=index){
			std::cout << "checkTopVoxels problem" << std::endl;
		}
		for(int j=0; j<10; j++){
			this->getRandomPositionAndOrientation(&pos, &angle, v, &rnd);
			if (this->getIndexOfTopLevelVoxel(pos)!=index){
				std::cout << "checkTopVoxels problem" << std::endl;
			}
		}

	}
}

Voxel *VoxelList::getVoxel(const RSAVector &pos, const RSAOrientation &angle){
	std::vector<Voxel *> *vTmp = this->voxelNeighbourGrid->getCell(pos);
	for(Voxel *v : *vTmp){
		if (v->isInside(pos, this->spatialVoxelSize, angle, this->angularVoxelSize)){
			return v;
		}
	}
	return nullptr;
}

Voxel *VoxelList::findVoxel(Voxel **list, size_t listSize, const RSAVector &pos, const RSAOrientation &angle){
	for(size_t i=0; i<listSize; i++){
		Voxel *v = list[i];
		if (v!=nullptr && v->isInside(pos, this->spatialVoxelSize, angle, this->angularVoxelSize)){
			return v;
		}
	}
	return nullptr;
}


void VoxelList::splitVoxel(Voxel *v, double spatialSize, double angularSize, Voxel **vRes){

	unsigned short spatialLoop = 1 << RSA_SPATIAL_DIMENSION;
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
		for(unsigned short j=0; j < RSA_SPATIAL_DIMENSION; j++){
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
		increment(inpos.data(), RSA_SPATIAL_DIMENSION, (unsigned char)1);
	} // for i
}


bool VoxelList::isVoxelInsidePacking(Voxel *v){
	RSAVector vpos = v->getPosition();
	RSAOrientation vangle = v->getOrientation();
	for(unsigned short i=0; i < RSA_SPATIAL_DIMENSION; i++){
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
		int arrayLenght = (int)round( pow(2, RSA_SPATIAL_DIMENSION+RSA_ANGULAR_DIMENSION) );
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


bool VoxelList::isTopLevelVoxelActive(Voxel *v){

	// checking if initial voxel containing v is active (do not have a shape inside)
	int index = this->getIndexOfTopLevelVoxel(v->getPosition());
	return this->activeTopLevelVoxels[index];
}


bool VoxelList::analyzeVoxel(Voxel *v, const RSAShape *s, RSABoundaryConditions *bc){

    if (!isTopLevelVoxelActive(v) || !this->isVoxelInsidePacking(v) )
		return true;

	return s->voxelInside(bc, v->getPosition(), v->getOrientation(), this->spatialVoxelSize, this->angularVoxelSize);
}


bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc, unsigned short depth){
	if (!this->disabled){ // && (depth > v->depth || depth==0) ){

	    if (!isTopLevelVoxelActive(v) || !this->isVoxelInsidePacking(v) )
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

		bool isInside = this->isVoxelInsideExclusionZone(v, this->spatialVoxelSize, this->angularVoxelSize, &shapes, bc, depth);

		v->depth = depth;
		v->lastAnalyzed = maxNo;
		return isInside;
	}
	return false;
}


size_t VoxelList::analyzeVoxels(RSABoundaryConditions *bc, NeighbourGrid<const RSAShape> *nl, unsigned short depth) {

	size_t begin = this->length;

	_OMP_PARALLEL_FOR
	for (size_t i = 0; i < this->length; i++) {
		Voxel *v = this->voxels[i];
		bool bRemove = this->analyzeVoxel(v, nl, bc, depth);
		if (bRemove){
			this->voxelNeighbourGrid->remove(v, v->getPosition());
			delete v;
			this->voxels[i] = nullptr;
		}
		if (i%10000 == 0){ std::cout << "."; std::cout.flush(); }
	}

	std::cout << " compacting" << std::flush;

	int endIndex = this->length-1;
	this->compactVoxelArray(this->voxels, endIndex);
	this->length = endIndex+1;

	return begin - this->length;
}


bool VoxelList::splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc){
	if (this->disabled)
		return false;
	size_t voxelsFactor = (size_t)round( pow(2, RSA_SPATIAL_DIMENSION+RSA_ANGULAR_DIMENSION) );
	if ((this->spatialVoxelSize<2*minDx && voxelsFactor*this->length > this->beginningVoxelNumber) || voxelsFactor*this->length > maxVoxels){
		return false;
	}

	int newListSize = voxelsFactor*(this->length);
	Voxel** newList = new Voxel*[ newListSize ];
	for(int i=0; i<newListSize; i++){
		newList[i] = nullptr;
	}

	Voxel ***aVoxels = new Voxel**[_OMP_MAXTHREADS];
	for(int i=0; i<_OMP_MAXTHREADS; i++){
		aVoxels[i] = new Voxel*[ voxelsFactor ];
	}

	_OMP_PARALLEL_FOR
	for(size_t i=0; i<this->length; i++){
		if (!this->analyzeVoxel(this->voxels[i], nl, bc)){ // dividing only not overlapping voxels
			this->splitVoxel(this->voxels[i], this->spatialVoxelSize/2.0, this->angularVoxelSize/2.0, aVoxels[_OMP_THREAD_ID]);
			for(size_t j=0; j<voxelsFactor; j++){
				Voxel *v = aVoxels[_OMP_THREAD_ID][j];
//				if(this->isVoxelInsidePacking(v) && ( nl==nullptr || bc==nullptr || !this->analyzeVoxel(v, nl, bc) ) ){
				if( nl==nullptr || bc==nullptr || !this->analyzeVoxel(v, nl, bc) ){
					if(this->voxels[i]->depth > 0){
						v->depth = this->voxels[i]->depth-1;
					}
					newList[i*voxelsFactor + j] = v;
				}else{
					delete aVoxels[_OMP_THREAD_ID][j];
				}
			}
		}
		delete this->voxels[i];
		if (i%10000 == 0){ std::cout << "." << std::flush; }
	}

	for(int i=0; i<_OMP_MAXTHREADS; i++){
		delete[] aVoxels[i];
	}
	delete[] aVoxels;

	delete[] this->voxels;

	int endIndex = newListSize - 1;

	std::cout << " compacting" << std::flush;

	this->compactVoxelArray(newList, endIndex);

	this->spatialVoxelSize = (this->spatialVoxelSize/2.0)*this->dxFactor;
	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);

	this->angularVoxelSize = (this->angularVoxelSize/2.0)*this->dxFactor;
	delete this->angularDistribution;
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);

	this->length = endIndex+1;
	this->voxels = newList;
	this->rebuildNeighbourGrid();

//	this->checkTopLevelVoxels();
	return true;
}


Voxel * VoxelList::getRandomVoxel(RND *rnd){
	double d = rnd->nextValue();
	return this->voxels[(int)(d*(this->length))];
}


Voxel * VoxelList::getVoxel(int i){
	return this->voxels[i];
}


void VoxelList::getRandomPositionAndOrientation(RSAVector *position, RSAOrientation *orientation, Voxel *v, RND *rnd){
	RSAVector vpos = v->getPosition();
	RSAOrientation vangle = v->getOrientation();

	for (ushort i=0; i < RSA_SPATIAL_DIMENSION; i++)
        (*position)[i] = vpos[i] + rnd->nextValue(this->spatialDistribution);
	for (ushort i=0; i < RSA_ANGULAR_DIMENSION; i++)
        (*orientation)[i] = vangle[i] + rnd->nextValue(this->angularDistribution);
}


double VoxelList::getVoxelSize(){
	return this->spatialVoxelSize;
}


double VoxelList::getVoxelAngularSize(){
	return this->angularVoxelSize;
}


Voxel* VoxelList::get(int i){
	return this->voxels[i];
}


size_t VoxelList::getLength() const{
	return this->length;
}


double VoxelList::getVoxelsSurface(){
	double result = 0;
	for(size_t i = 0; i< this->length; i++){
		double s = 1.0;
		Voxel *v = this->voxels[i];
		RSAVector position = v->getPosition();
		for(unsigned short j=0; j<RSA_SPATIAL_DIMENSION; j++){
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


std::string VoxelList::toPovray(){
	std::string sRes = "";

	for(size_t i=0; i<this->length; i++){
		sRes += this->voxels[i]->toPovray(this->spatialVoxelSize) + "\r\n";
	}
	return sRes;
}


std::string VoxelList::toWolfram(){
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

	this->voxelNeighbourGrid->clear();
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
	this->rebuildNeighbourGrid();
}
