/*
 * VoxelList.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#include <iostream>
#include <cmath>
#include <algorithm>

#include "Voxel.h"
#include "math.h"
#include "NeighbourGrid.h"
#include "Shape.h"

#include "Utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif


/**
 * d - requested initial size of a voxel
 */
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::findInitialVoxelSize(double d){
	double dRes = 1.0;
	while (dRes > d)
		dRes /= 2;
	while (2*dRes < d)
		dRes *= 2;
	return dRes;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::findInitialVoxelAngularSize(double d){
	double dRes = 4.0;
	while (dRes < d)
		dRes *= 2;
	while (dRes > d)
		dRes /= 2;
	return 2*dRes;
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
inline int VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getLinearNumberOfVoxels(double vs){
	return (int)(this->size/vs) + 1;
}

/**
dim - packing dimension
s - packing size (linear)
d - requested initial size of a voxel
**/
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::VoxelList(double s, double d, double ad){
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
	int voxelsLength = (int)(pow(n, SPATIAL_DIMENSION)+0.5);
	this->voxels = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>*[voxelsLength];
	this->activeTopLevelVoxels = new bool[voxelsLength];
	this->voxelNeighbourGrid = new NeighbourGrid<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>>(SPATIAL_DIMENSION, this->voxelSize*n, n);


	int in[SPATIAL_DIMENSION];
	for(ushort i=0; i<SPATIAL_DIMENSION; i++)
		in[i] = 0;
	int index = 0;
	do{
		std::copy(in, in+SPATIAL_DIMENSION, offset[index]);
		index++;
	}while(increment(in, SPATIAL_DIMENSION, 1));

	this->initVoxels();
	this->voxelSize *= this->dxFactor;
	this->beginningVoxelNumber = voxelsLength;
	this->last = voxelsLength-1;
	this->fillNeighbourGrid();

//		this.checkIndexes();
	}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::~VoxelList() {
	delete this->spatialDistribution;

	for(int i=0; i<=this->last; i++){
		delete this->voxels[i];
	}
	delete[] this->voxels;
	delete[] this->activeTopLevelVoxels;
	delete this->voxelNeighbourGrid;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::disable(){
	this->disabled = true;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::initVoxels(){
	int n = this->getLinearNumberOfVoxels(this->voxelSize);
	double position[SPATIAL_DIMENSION];
	double orientation[ANGULAR_DIMENSION];
	int in[SPATIAL_DIMENSION];

	for(ushort i = 0; i<SPATIAL_DIMENSION; i++)
		in[i] = 0;

	for(ushort i = 0; i<ANGULAR_DIMENSION; i++)
		orientation[i] = 0;


	int i, index = 0;
	do{
		for(unsigned char i=0; i<SPATIAL_DIMENSION; i++){
			position[i] = this->voxelSize*in[i]; // da point to the "left bottom" corner of a voxel
		}
		i = position2i(position, SPATIAL_DIMENSION, n*this->voxelSize, this->voxelSize, n);
		if(index!=i){
			std::cout << "VoxelList::initVoxels: Problem: " << index << " != " << i << std::endl;
		}

		this->voxels[index] = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>(position, orientation);
		this->voxels[index]->index = index;
		this->activeTopLevelVoxels[index] = true;
		index++;
	}while(increment(in, SPATIAL_DIMENSION, n-1));
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::fillNeighbourGrid(){
	this->voxelNeighbourGrid->clear();
	for(int i=0; i<=this->last; i++){
		this->voxelNeighbourGrid->add(this->voxels[i], this->voxels[i]->getPosition());
	}
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getNeighbours(std::vector<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *result, double *da){
	return this->voxelNeighbourGrid->getNeighbours(result, da);
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::checkIndexes(){
	for(int i=0; i<=this->last; i++){
		if(this->voxels[i]->index!=i)
			std::cout << "VoxelList::checkIndexes: Error " << i << std::endl;
	}
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::remove(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v){
	if (this->disabled)
		return;

	int index = v->index;

	if (index!=last){
		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *vl = this->voxels[last];
		vl->index = index;
		this->voxels[index] = vl;
	}
	this->voxels[last] = NULL;
	this->last--;
	this->voxelNeighbourGrid->remove(v, v->getPosition());
	delete v;

//		this.checkIndexes();
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::compactVoxelArray(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **list, int &endIndex){
	int beginIndex = 0;

	while(list[endIndex]==NULL && endIndex > -1)
		(endIndex)--;

	while (beginIndex<endIndex){
		if(list[beginIndex]==NULL){
			list[beginIndex] = list[endIndex];
			list[beginIndex]->index = beginIndex;
			list[endIndex] = NULL;
		}
		beginIndex++;
		while(list[endIndex]==NULL && endIndex > -1)
			endIndex--;
	}
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::remove(std::vector<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *vVoxels){
	if (this->disabled)
		return;
	for(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v : *vVoxels){
		this->voxels[v->index] = NULL;
		delete v;
	}

	int endIndex = this->last;
	this->compactVoxelArray(this->voxels, endIndex);

	this->voxelNeighbourGrid->clear();
	this->last = endIndex;
	this->fillNeighbourGrid();
//	this->checkIndexes();

}



template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getIndexOfTopLevelVoxel(double *da){

	int n = (int)(this->size/this->initialVoxelSize) + 1;
	int index = position2i(da, SPATIAL_DIMENSION, n*this->initialVoxelSize, this->initialVoxelSize, n);
	return index;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::removeTopLevelVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v){
	if (this->disabled)
		return;
	this->activeTopLevelVoxels[this->getIndexOfTopLevelVoxel(v->getPosition())]=false;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::checkTopLevelVoxels(){

	RND rnd;
	double pos[SPATIAL_DIMENSION];
	double angle[ANGULAR_DIMENSION];

	for(int i=0; i<=this->last; i++){
		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v = this->voxels[i];
		this->getRandomPositionAndOrientation(pos, angle, v, &rnd);
		int index = this->getIndexOfTopLevelVoxel(pos);
		if (this->getIndexOfTopLevelVoxel(v->getPosition())!=index){
			std::cout << "checkTopVoxels problem" << std::endl;
		}
		for(int j=0; j<10; j++){
			this->getRandomPositionAndOrientation(pos, angle, v, &rnd);
			if (this->getIndexOfTopLevelVoxel(pos)!=index){
				std::cout << "checkTopVoxels problem" << std::endl;
			}
		}

	}
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> * VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVoxel(double *pos, const double *angle){
	std::vector<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *vTmp = this->voxelNeighbourGrid->getCell(pos);
	for(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v : *vTmp){
		if (v->isInside(pos, this->voxelSize, angle, this->angularVoxelSize)){
			return v;
		}
	}
	return NULL;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::splitVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, double spatialSize, double angularSize, Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **vRes){

	unsigned short spatialLoop = 1 << SPATIAL_DIMENSION;
	unsigned short angularLoop = 1 << ANGULAR_DIMENSION;

	int inpos[SPATIAL_DIMENSION];
	double position[SPATIAL_DIMENSION];
	int inangle[ANGULAR_DIMENSION];
	double orientation[ANGULAR_DIMENSION];


	double* vpos = v->getPosition();
	for(ushort j=0; j < SPATIAL_DIMENSION; j++){
		inpos[j] = 0;
		position[j] = vpos[j];
	}

	double *vangle = v->getOrientation();
	for(ushort j=0; j < ANGULAR_DIMENSION; j++){
		inangle[j] = 0;
		orientation[j] = vangle[j];
	}

	for(unsigned short i=0; i<spatialLoop; i++){
		for(unsigned short j=0; j < SPATIAL_DIMENSION; j++){
			position[j] = vpos[j] + inpos[j]*spatialSize;
		}
		for(unsigned short j=0; j < ANGULAR_DIMENSION; j++){
			inangle[j] = 0;
		}
		for(unsigned short j=0; j<angularLoop; j++){
			for(unsigned short k=0; k<ANGULAR_DIMENSION; k++){
				orientation[k] = vangle[k] + inangle[k]*angularSize;
			}
			vRes[i*angularLoop + j] = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>(position, orientation);
			increment(inangle, ANGULAR_DIMENSION, (unsigned char)1);
		} // for j
		increment(inpos, SPATIAL_DIMENSION, (unsigned char)1);
	} // for i
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::isVoxelInsidePacking(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v){
	double* vpos = v->getPosition();
	double *vangle = v->getOrientation();
	for(unsigned short i=0; i < SPATIAL_DIMENSION; i++){
		if (vpos[i] >= this->size){
			return false;
		}
	}
	for(unsigned short i=0; i<ANGULAR_DIMENSION; i++){
		if (vangle[i] >= this->angularSize){
			return false;
		}
	}
	return true;
}

/**
 * returns true when the whole voxel is inside an exclusion area of the shape s
 */
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::isVoxelInsideExclusionZone(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, double spatialSize, double angularSize, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s, BoundaryConditions *bc){

	double* vpos = v->getPosition();
	double position[SPATIAL_DIMENSION];
	int counterSize = 1 << SPATIAL_DIMENSION;
	bool shapeCoversVertices = true;
	for(int i=0; i<counterSize; i++){
		for(ushort j=0; j<SPATIAL_DIMENSION; j++){
			position[j] = vpos[j] + this->offset[i][j]*spatialSize;
		}
		if(!s->pointInside(bc, position, v->getOrientation(), angularSize)){
			shapeCoversVertices = false;
			break;
		}
	}
	return shapeCoversVertices;
}

/**
 * returns true when the whole voxel is inside an exclusion area of any shape in shapes
 * To determine it the method tires to split voxel up to level of maxDepth
 */
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::isVoxelInsideExclusionZone(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, double spatialSize, double angularSize, std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *shapes, BoundaryConditions *bc, unsigned short depth){
	// if voxel is outside the packing it is inside exclusion zone
	if (!this->isVoxelInsidePacking(v))
		return true;
	// otherwise checking

	bool isInside = false;
	for(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s : *shapes){
		isInside = this->isVoxelInsideExclusionZone(v, spatialSize, angularSize, s, bc);
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
		int arrayLenght = (int)round( pow(2, SPATIAL_DIMENSION+ANGULAR_DIMENSION) );
		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **aVoxels = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>*[ arrayLenght ];
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

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s, BoundaryConditions *bc){

	double *vpos = v->getPosition();
	// checking if initial voxel containing v is active (do not have a shape inside)
	int index = this->getIndexOfTopLevelVoxel(vpos);
	if(this->activeTopLevelVoxels[index]==false)
		return true;


	return this->isVoxelInsideExclusionZone(v, this->voxelSize, this->angularVoxelSize, s, bc);
}

/*
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *neighbours, BoundaryConditions *bc){
	double* vpos = v->getPosition();
	// checking if initial voxel containing v is active (do not have a shape inside)
	int index = this->getIndexOfTopLevelVoxel(vpos);
	if(this->activeTopLevelVoxels[index]==false)
		return true;

	std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *vNeighbours;

	return this->isVoxelInsideExclusionZone(v, this->voxelSize, this->angularVoxelSize, neighbours, bc);
}
*/

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc, unsigned short depth){
	if (!this->disabled && (depth > v->depth || depth==0) ){

		double* vpos = v->getPosition();
		// checking if initial voxel containing v is active (do not have a shape inside)
		int index = this->getIndexOfTopLevelVoxel(vpos);
		if(this->activeTopLevelVoxels[index]==false)
			return true;

		std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> tmpShapes, shapes;
		nl->getNeighbours(&tmpShapes, v->getPosition());

		int maxNo = v->lastAnalyzed;
		for(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s: tmpShapes){
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

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::splitVoxels(double minDx, int maxVoxels, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc){
	if (this->disabled)
		return false;
	int voxelsFactor = (int)round( pow(2, SPATIAL_DIMENSION+ANGULAR_DIMENSION) );
	if ((this->voxelSize<2*minDx && voxelsFactor*this->last > this->beginningVoxelNumber) || voxelsFactor*this->last > maxVoxels){
		return false;
	}

	int newListSize = voxelsFactor*(this->last+1);
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>** newList = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>*[ newListSize ];
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
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> ***aVoxels = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>**[maxthreads];
	for(int i=0; i<maxthreads; i++){
		aVoxels[i] = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>*[ voxelsFactor ];
	}

	#pragma omp parallel for
	for(int i=0; i<=this->last; i++){
		int tid = omp_get_thread_num();
		this->splitVoxel(this->voxels[i], this->voxelSize, this->angularVoxelSize, aVoxels[tid]);
		for(int j=0; j<voxelsFactor; j++){
			Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v = aVoxels[tid][j];
			if(this->isVoxelInsidePacking(v) && ( nl==NULL || bc==NULL || !this->analyzeVoxel(v, nl, bc) ) ){
					v->index = i*voxelsFactor + j;
					newList[i*voxelsFactor + j] = v;
			}else{
				delete aVoxels[tid][j];
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

	this->last = endIndex;
	this->voxels = newList;
	this->fillNeighbourGrid();
//	this->checkIndexes();
//	this->checkTopLevelVoxels();
	return true;
}

#else

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::splitVoxels(double minDx, int maxVoxels, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc){
	if (this->disabled)
		return false;
	int voxelsFactor = (int)round( pow(2, SPATIAL_DIMENSION+ANGULAR_DIMENSION) );
	if ((this->voxelSize<2*minDx && voxelsFactor*this->last > this->beginningVoxelNumber) || voxelsFactor*this->last > maxVoxels){
		return false;
	}

	int newListSize = voxelsFactor*(this->last+1);
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>** newList = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>*[ newListSize ];

	this->voxelSize = (this->voxelSize/2.0)*this->dxFactor;
	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);

	this->angularVoxelSize = this->angularVoxelSize/2.0;
	delete this->angularDistribution;
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);

	int index = 0;

	int inpos[SPATIAL_DIMENSION];
	double position[SPATIAL_DIMENSION];
	int inangle[ANGULAR_DIMENSION];
	double orientation[ANGULAR_DIMENSION];
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **aVoxels = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>*[voxelsFactor];

	for(int i=0; i<=this->last; i++){
		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v = this->voxels[i];

		this->splitVoxel(this->voxels[i], this->voxelSize, this->angularVoxelSize, aVoxels);
		for(int j=0; j<voxelsFactor; j++){
			if(this->isVoxelInsidePacking(aVoxels[j]) && ( nl==NULL || bc==NULL || !this->analyzeVoxel(aVoxels[j], nl, bc) ) ){
					aVoxels[j]->index = index;
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

	this->last = index-1;
	this->voxels = newList;
	this->fillNeighbourGrid();
//	this->checkIndexes();
//	this->checkTopLevelVoxels();
	return true;
}

#endif


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> * VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getRandomVoxel(RND *rnd){
	double d = rnd->nextValue();
	return this->voxels[(int)(d*(this->last+1))];
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> * VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVoxel(int i){
	return this->voxels[i];
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getRandomPositionAndOrientation(double *position, double *orientation, Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, RND *rnd){
	double *vpos = v->getPosition();
	double *vangle = v->getOrientation();

	for (ushort i=0; i < SPATIAL_DIMENSION; i++)
		position[i] = vpos[i] + rnd->nextValue(this->spatialDistribution);
	for (ushort i=0; i < ANGULAR_DIMENSION; i++)
		orientation[i] = vangle[i] + rnd->nextValue(this->angularDistribution);
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVoxelSize(){
	return this->voxelSize;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVoxelAngularSize(){
	return this->angularVoxelSize;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>* VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::get(int i){
	return this->voxels[i];
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::length() const{
	return this->last+1;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVoxelsSurface(){
	return (this->last+1)*pow(this->voxelSize, SPATIAL_DIMENSION)*pow((this->angularVoxelSize/(this->initialAngularVoxelSize)), ANGULAR_DIMENSION);
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toPovray(){
	std::string sRes = "";

	for(int i=0; i<=this->last; i++){
		sRes += this->voxels[i]->toPovray(this->voxelSize) + "\r\n";
	}
	return sRes;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toWolfram(){
	std::stringstream out;

	for(int i=0; i<=this->last; i++){
		out << this->voxels[i]->toWolfram(this->voxelSize, this->angularVoxelSize);
		if (i!=this->last)
			out << ", ";
		out << std::endl;
		if (ANGULAR_DIMENSION > 0){
			out << "(* angles: [ " << this->voxels[i]->getOrientation()[0] << ", " << (this->voxels[i]->getOrientation()[0] + this->angularVoxelSize) << ") *)" << std::endl;
		}
	}
	return out.str();
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::store(std::ostream &f) const{
	unsigned short sd = SPATIAL_DIMENSION;
	unsigned short ad = ANGULAR_DIMENSION;
	f.write((char *)(&sd), sizeof(unsigned char));
	if (ad>0)
		f.write((char *)(&ad), sizeof(unsigned char));
	f.write((char *) &this->voxelSize, sizeof(double));
	f.write((char *) &this->angularVoxelSize, sizeof(double));

	int size = this->length();
	f.write((char *)(&size), sizeof(int));
	for(int i=0; i<size; i++){
		this->voxels[i]->store(f);
	}
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::restore(std::istream &f){
	unsigned char sd = SPATIAL_DIMENSION;
	unsigned char ad = ANGULAR_DIMENSION;

	f.read((char *)(&sd), sizeof(unsigned char));
	if (ad > 0)
		f.read((char *)(&ad), sizeof(unsigned char));

	if (sd!=SPATIAL_DIMENSION || ad!=ANGULAR_DIMENSION){
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
	for(int i=0; i<=this->last; i++){
		delete this->voxels[i];
	}
	delete[] this->voxels;

	int size;
	int topIndex;
	f.read((char *)&size, sizeof(int));
	this->voxels = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[size];
	for(int i=0; i<size; i++){
		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>();
		v->restore(f);
		this->voxels[i] = v;
		v->index = i;
		topIndex = this->getIndexOfTopLevelVoxel(v->getPosition());
		this->activeTopLevelVoxels[topIndex] = true;
	}
	this->last = size-1;
	this->fillNeighbourGrid();
}
