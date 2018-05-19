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

double VoxelList::findInitialVoxelSize(double d){
	double dRes = 1.0;
	while (dRes > d)
		dRes /= 2;
	while (2*dRes < d)
		dRes *= 2;
	return dRes;
}


double VoxelList::findInitialVoxelAngularSize(double d){
	double dRes = 4.0;
	while (dRes < d)
		dRes *= 2;
	while (dRes > d)
		dRes /= 2;
	return 2*dRes;
}



inline int VoxelList::getLinearNumberOfVoxels(double vs){
	return (int)(this->size/vs) + 1;
}

/**
dim - packing dimension
s - packing size (linear)
d - requested initial size of a voxel
**/

VoxelList::VoxelList(double s, double d, double ad){
	this->size = s;
	this->voxelSize = this->findInitialVoxelSize(d);
	this->initialVoxelSize = this->voxelSize;

	this->angularSize = ad;
	this->angularVoxelSize = this->findInitialVoxelAngularSize(ad);
	this->initialAngularVoxelSize = this->angularVoxelSize;

	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);
	this->disabled = false;

	int n = this->getLinearNumberOfVoxels(this->voxelSize);
	int voxelsLength = (int)(pow(n, RSA_SPATIAL_DIMENSION)+0.5);
	this->voxels = new Voxel*[voxelsLength];
	this->activeTopLevelVoxels = new bool[voxelsLength];
	this->voxelNeighbourGrid = new NeighbourGrid<Voxel>(RSA_SPATIAL_DIMENSION, this->voxelSize*n, n);

	this->initVoxels();
	this->voxelSize *= this->dxFactor;
	this->beginningVoxelNumber = voxelsLength;
	this->length = voxelsLength;
	this->fillNeighbourGrid();

//		this.checkIndexes();
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
	int n = this->getLinearNumberOfVoxels(this->voxelSize);
	double position[RSA_SPATIAL_DIMENSION];
	std::array<double, RSA_ANGULAR_DIMENSION> orientation;
	int in[RSA_SPATIAL_DIMENSION];

	for(ushort i = 0; i<RSA_SPATIAL_DIMENSION; i++)
		in[i] = 0;

	for(ushort i = 0; i<RSA_ANGULAR_DIMENSION; i++)
		orientation[i] = 0;


	int index = 0;
	do{
		for(unsigned char i=0; i<RSA_SPATIAL_DIMENSION; i++){
			position[i] = this->voxelSize*in[i]; // da point to the "left bottom" corner of a voxel
		}
/*
		int i = position2i(position, RSA_SPATIAL_DIMENSION, n*this->voxelSize, this->voxelSize, n);
		if(index!=i){
			std::cout << "VoxelList::initVoxels: Problem: " << index << " != " << i << std::endl;
		}
*/
		this->voxels[index] = new Voxel(position, orientation);
		this->activeTopLevelVoxels[index] = true;
		index++;
	}while(increment(in, RSA_SPATIAL_DIMENSION, n-1));
}


void VoxelList::fillNeighbourGrid(){
	this->voxelNeighbourGrid->clear();
	for(size_t i=0; i<this->length; i++){
		this->voxelNeighbourGrid->add(this->voxels[i], this->voxels[i]->getPosition());
	}
}


void VoxelList::getNeighbours(std::vector<Voxel *> *result, double *da){
	return this->voxelNeighbourGrid->getNeighbours(result, da);
}

void VoxelList::compactVoxelArray(Voxel **list, int &endIndex){
	int beginIndex = 0;

	while(list[endIndex]==NULL && endIndex > -1)
		(endIndex)--;

	while (beginIndex<endIndex){
		if(list[beginIndex]==NULL){
			list[beginIndex] = list[endIndex];
			list[endIndex] = NULL;
		}
		beginIndex++;
		while(list[endIndex]==NULL && endIndex > -1)
			endIndex--;
	}
}

int VoxelList::getIndexOfTopLevelVoxel(double *da){

	int n = (int)(this->size/this->initialVoxelSize) + 1;
	int index = position2i(da, RSA_SPATIAL_DIMENSION, n*this->initialVoxelSize, this->initialVoxelSize, n);
	return index;
}


void VoxelList::removeTopLevelVoxel(Voxel *v){
	if (this->disabled)
		return;
	this->activeTopLevelVoxels[this->getIndexOfTopLevelVoxel(v->getPosition())]=false;
}


void VoxelList::checkTopLevelVoxels(){

	RND rnd;
	double pos[RSA_SPATIAL_DIMENSION];
	std::array<double, RSA_ANGULAR_DIMENSION> angle;

	for(size_t i=0; i<this->length; i++){
		Voxel *v = this->voxels[i];
		this->getRandomPositionAndOrientation(pos, angle.data(), v, &rnd);
		int index = this->getIndexOfTopLevelVoxel(pos);
		if (this->getIndexOfTopLevelVoxel(v->getPosition())!=index){
			std::cout << "checkTopVoxels problem" << std::endl;
		}
		for(int j=0; j<10; j++){
			this->getRandomPositionAndOrientation(pos, angle.data(), v, &rnd);
			if (this->getIndexOfTopLevelVoxel(pos)!=index){
				std::cout << "checkTopVoxels problem" << std::endl;
			}
		}

	}
}


Voxel * VoxelList::getVoxel(double *pos, const std::array<double, RSA_ANGULAR_DIMENSION> &angle){
	std::vector<Voxel *> *vTmp = this->voxelNeighbourGrid->getCell(pos);
	for(Voxel *v : *vTmp){
		if (v->isInside(pos, this->voxelSize, angle, this->angularVoxelSize)){
			return v;
		}
	}
	return NULL;
}


void VoxelList::splitVoxel(Voxel *v, double spatialSize, double angularSize, Voxel **vRes){

	unsigned short spatialLoop = 1 << RSA_SPATIAL_DIMENSION;
	unsigned short angularLoop = 1 << RSA_ANGULAR_DIMENSION;

	int inpos[RSA_SPATIAL_DIMENSION];
	double position[RSA_SPATIAL_DIMENSION];
	std::array<int, RSA_ANGULAR_DIMENSION> inangle;
	std::array<double, RSA_ANGULAR_DIMENSION> orientation;


	double* vpos = v->getPosition();
	for(ushort j=0; j < RSA_SPATIAL_DIMENSION; j++){
		inpos[j] = 0;
		position[j] = vpos[j];
	}

	std::array<double, RSA_ANGULAR_DIMENSION> vangle = v->getOrientation();
    orientation = vangle;
	for(ushort j=0; j < RSA_ANGULAR_DIMENSION; j++)
		inangle[j] = 0;

	for(unsigned short i=0; i<spatialLoop; i++){
		for(unsigned short j=0; j < RSA_SPATIAL_DIMENSION; j++){
			position[j] = vpos[j] + inpos[j]*spatialSize;
		}
		for(unsigned short j=0; j < RSA_ANGULAR_DIMENSION; j++){
			inangle[j] = 0;
		}
		for(unsigned short j=0; j<angularLoop; j++){
			for(unsigned short k=0; k<RSA_ANGULAR_DIMENSION; k++){
				orientation[k] = vangle[k] + inangle[k]*angularSize;
			}
			vRes[i*angularLoop + j] = new Voxel(position, orientation);
			increment(inangle.data(), RSA_ANGULAR_DIMENSION, (unsigned char)1);
		} // for j
		increment(inpos, RSA_SPATIAL_DIMENSION, (unsigned char)1);
	} // for i
}


bool VoxelList::isVoxelInsidePacking(Voxel *v){
	double* vpos = v->getPosition();
	std::array<double, RSA_ANGULAR_DIMENSION> vangle = v->getOrientation();
	for(unsigned short i=0; i < RSA_SPATIAL_DIMENSION; i++){
		if (vpos[i] >= this->size){
			return false;
		}
	}
	for(unsigned short i=0; i<RSA_ANGULAR_DIMENSION; i++){
		if (vangle[i] >= this->angularSize){
			return false;
		}
	}
	return true;
}

/**
 * returns true when the whole voxel is inside an exclusion area of any shape in shapes
 * To determine it the method tires to split voxel up to level of maxDepth
 */

bool VoxelList::isVoxelInsideExclusionZone(Voxel *v, double spatialSize, double angularSize, std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *shapes, BoundaryConditions *bc, unsigned short depth){
	// if voxel is outside the packing it is inside exclusion zone
	if (!this->isVoxelInsidePacking(v))
		return true;
	// otherwise checking

	bool isInside = false;
	for(Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s : *shapes){
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

	double *vpos = v->getPosition();
	// checking if initial voxel containing v is active (do not have a shape inside)
	int index = this->getIndexOfTopLevelVoxel(vpos);
	return this->activeTopLevelVoxels[index];
}


bool VoxelList::analyzeVoxel(Voxel *v, Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s, BoundaryConditions *bc){

	if (!isTopLevelVoxelActive(v))
		return true;

	return s->voxelInside(bc, v->getPosition(), v->getOrientation(), this->voxelSize, this->angularVoxelSize);
}


bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc, unsigned short depth){
	if (!this->disabled && (depth > v->depth || depth==0) ){

	if (!isTopLevelVoxelActive(v))
		return true;

	std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> tmpShapes, shapes;
		nl->getNeighbours(&tmpShapes, v->getPosition());

		int maxNo = v->lastAnalyzed;
		for(Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s: tmpShapes){
			if (v->lastAnalyzed < s->no || depth > v->depth){
				shapes.push_back(s);
				if(maxNo < s->no)
					maxNo = s->no;
			}
		}

		bool isInside = this->isVoxelInsideExclusionZone(v, this->voxelSize, this->angularVoxelSize, &shapes, bc, depth);

		v->depth = depth;
		v->lastAnalyzed = maxNo;
		return isInside;
	}
	return false;
}

#ifdef _OPENMP

size_t VoxelList::analyzeVoxels(BoundaryConditions *bc, NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> *nl, unsigned short depth) {

	size_t begin = this->length;

	size_t newListSize = this->length;
	Voxel** newList = new Voxel*[ newListSize ];

	#pragma omp parallel for
	for(size_t i=0; i<newListSize; i++){
		newList[i] = NULL;
	}

	size_t newIndex = 0;
	#pragma omp parallel for
	for (size_t i = 0; i < this->length; i++) {
		Voxel *v = this->voxels[i];
		bool bRemove = this->analyzeVoxel(v, nl, bc, depth);
		if (!bRemove) {
			#pragma omp critical
			{
			  newList[newIndex] = v;
			  newIndex++;
		    }
		}else{
			delete v;
		}
		if (i%10000 == 0){ std::cout << "."; std::cout.flush(); }
	}

	delete[] this->voxels;

	this->length = newIndex;
	this->voxels = newList;

	return begin - this->length;
}

#else

size_t VoxelList::analyzeVoxels(BoundaryConditions *bc, NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> *nl, unsigned short depth) {

	size_t begin = this->length;

	size_t newListSize = this->length;
	Voxel** newList = new Voxel*[ newListSize ];
	for(size_t i=0; i<newListSize; i++){
		newList[i] = NULL;
	}

	size_t newIndex = 0;
	for (size_t i = 0; i < this->length; i++) {
		Voxel *v = this->voxels[i];
		bool bRemove = this->analyzeVoxel(v, nl, bc, depth);
		if (!bRemove) {
			newList[newIndex] = v;
			newIndex++;
		}else{
			delete v;
		}
		if (i%10000 == 0){ std::cout << "."; std::cout.flush(); }
	}

	delete[] this->voxels;

	this->length = newIndex;
	this->voxels = newList;

	return begin - this->length;
}

#endif


#ifdef _OPENMP


bool VoxelList::splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc){
	if (this->disabled)
		return false;
	size_t voxelsFactor = (size_t)round( pow(2, RSA_SPATIAL_DIMENSION+RSA_ANGULAR_DIMENSION) );
	if ((this->voxelSize<2*minDx && voxelsFactor*this->length > this->beginningVoxelNumber) || voxelsFactor*this->length > maxVoxels){
		return false;
	}

	int newListSize = voxelsFactor*(this->length);
	Voxel** newList = new Voxel*[ newListSize ];
	for(int i=0; i<newListSize; i++){
		newList[i] = NULL;
	}

	this->voxelSize = (this->voxelSize/2.0)*this->dxFactor;
	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);

	this->angularVoxelSize = (this->angularVoxelSize/2.0)*this->dxFactor;
	delete this->angularDistribution;
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);

	int maxthreads = omp_get_max_threads();
	Voxel ***aVoxels = new Voxel**[maxthreads];
	for(int i=0; i<maxthreads; i++){
		aVoxels[i] = new Voxel*[ voxelsFactor ];
	}

	#pragma omp parallel for
	for(size_t i=0; i<this->length; i++){
		int tid = omp_get_thread_num();
		if (!this->analyzeVoxel(this->voxels[i], nl, bc)){ // dividing only not overlapping voxels
			this->splitVoxel(this->voxels[i], this->voxelSize, this->angularVoxelSize, aVoxels[tid]);
			for(size_t j=0; j<voxelsFactor; j++){
				Voxel *v = aVoxels[tid][j];
				if(this->isVoxelInsidePacking(v) && ( nl==NULL || bc==NULL || !this->analyzeVoxel(v, nl, bc) ) ){
					if(this->voxels[i]->depth > 0){
						v->depth = this->voxels[i]->depth-1;
					}
					newList[i*voxelsFactor + j] = v;
				}else{
					delete aVoxels[tid][j];
				}
			}
		}
		delete this->voxels[i];
		if (i%10000 == 0){ std::cout << "." << std::flush; }
	}

	for(int i=0; i<maxthreads; i++){
		delete[] aVoxels[i];
	}
	delete[] aVoxels;

	delete[] this->voxels;

	int endIndex = newListSize - 1;;

	std::cout << " compacting" << std::flush;

	this->compactVoxelArray(newList, endIndex);

	this->voxelNeighbourGrid->clear();

	this->length = endIndex+1;
	this->voxels = newList;
	this->fillNeighbourGrid();
//	this->checkIndexes();
//	this->checkTopLevelVoxels();
	return true;
}

#else


bool VoxelList::splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc){
	if (this->disabled)
		return false;
	int voxelsFactor = (int)round( pow(2, RSA_SPATIAL_DIMENSION+RSA_ANGULAR_DIMENSION) );
	if ((this->voxelSize<2*minDx && voxelsFactor*this->length > this->beginningVoxelNumber) || voxelsFactor*this->length > maxVoxels){
		return false;
	}

	int newListSize = voxelsFactor*(this->length);
	Voxel** newList = new Voxel*[ newListSize ];

	this->voxelSize = (this->voxelSize/2.0)*this->dxFactor;
	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);

	this->angularVoxelSize = this->angularVoxelSize/2.0;
	delete this->angularDistribution;
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);

	int index = 0;

	Voxel **aVoxels = new Voxel*[voxelsFactor];

	for(size_t i=0; i<this->length; i++){

		this->splitVoxel(this->voxels[i], this->voxelSize, this->angularVoxelSize, aVoxels);
		for(int j=0; j<voxelsFactor; j++){
			if(this->isVoxelInsidePacking(aVoxels[j]) && ( nl==NULL || bc==NULL || !this->analyzeVoxel(aVoxels[j], nl, bc) ) ){
					if(this->voxels[i]->depth > 0){
						aVoxels[j]->depth = this->voxels[i]->depth-1;
					}
					newList[index] = aVoxels[j];
					index++;
			}else{
				delete aVoxels[j];
			}
		}
		delete this->voxels[i];
		if (i%10000 == 0){ std::cout << "."; std::cout.flush(); }
	}

	delete[] this->voxels;
	this->voxelNeighbourGrid->clear();

	this->length = index;
	this->voxels = newList;
	this->fillNeighbourGrid();
//	this->checkIndexes();
//	this->checkTopLevelVoxels();
	return true;
}

#endif



Voxel * VoxelList::getRandomVoxel(RND *rnd){
	double d = rnd->nextValue();
	return this->voxels[(int)(d*(this->length))];
}


Voxel * VoxelList::getVoxel(int i){
	return this->voxels[i];
}


void VoxelList::getRandomPositionAndOrientation(double *position, double *orientation, Voxel *v, RND *rnd){
	double *vpos = v->getPosition();
	std::array<double, RSA_ANGULAR_DIMENSION> vangle = v->getOrientation();

	for (ushort i=0; i < RSA_SPATIAL_DIMENSION; i++)
		position[i] = vpos[i] + rnd->nextValue(this->spatialDistribution);
	for (ushort i=0; i < RSA_ANGULAR_DIMENSION; i++)
		orientation[i] = vangle[i] + rnd->nextValue(this->angularDistribution);
}


double VoxelList::getVoxelSize(){
	return this->voxelSize;
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
	return (this->length)*pow(this->voxelSize, RSA_SPATIAL_DIMENSION)*pow((this->angularVoxelSize/(this->initialAngularVoxelSize)), RSA_ANGULAR_DIMENSION);
}


std::string VoxelList::toPovray(){
	std::string sRes = "";

	for(size_t i=0; i<this->length; i++){
		sRes += this->voxels[i]->toPovray(this->voxelSize) + "\r\n";
	}
	return sRes;
}


std::string VoxelList::toWolfram(){
	std::stringstream out;

	for(size_t i=0; i<this->length; i++){
		out << this->voxels[i]->toWolfram(this->voxelSize, this->angularVoxelSize);
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
	f.write((char *) &this->voxelSize, sizeof(double));
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
	f.read((char *) &(this->voxelSize), sizeof(double));
	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);

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
	this->fillNeighbourGrid();
}