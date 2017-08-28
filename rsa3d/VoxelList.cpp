/*
 * VoxelList.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */


#include "VoxelList.h"

#include <iostream>
#include <cmath>
#include <algorithm>

#include "Voxel.h"
#include "math.h"
#include "NeighbourGrid.h"
#include "Shape.h"

#include "Utils.h"

/**
dim - packing dimension
s - packing size (linear)
d - requested initial size of a voxel
**/
VoxelList::VoxelList(unsigned char dim, double s, double d){
	this->dimension = dim;
	this->size = s;
	int n = (int)(s/d)+1;
	this->voxelSize = (s/n); // actual initial size of a voxel (should be not greater than d)
	this->initialVoxelSize = this->voxelSize;

	this->distribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);
	this->disabled = false;

	n = (int)(this->size/this->voxelSize + 0.5);
	int voxelsLength = (int)(pow(n, dim)+0.5);
	this->voxels = new Voxel*[voxelsLength];
	this->activeTopLevelVoxels = new bool[voxelsLength];
	this->voxelNeighbourGrid = new NeighbourGrid(dim, s, this->voxelSize);

	counter = new double*[(1 << this->dimension)]; // matrix of d-dimensional offsets to 2^d voxel vertices

	int *in = new int[this->dimension];
	for(unsigned char i=0; i<this->dimension; i++){
		in[i] = 0;
	}
	int index = 0;
	do{
		counter[index] = new double[this->dimension];

		for(unsigned char i=0; i<this->dimension; i++){
			counter[index][i] = in[i];
		}
		index++;
	}while(increment(in, this->dimension, 1));
	delete[] in;


	this->initVoxels(dim);
	this->voxelSize *= this->dxFactor;
	this->beginningVoxelNumber = voxelsLength;
	this->last = voxelsLength-1;
	this->fillNeighbourGrid();

//		this.checkIndexes();
	}

VoxelList::~VoxelList() {
	delete this->distribution;

	for(int i=0; i<=this->last; i++){
		delete this->voxels[i];
	}
	delete[] this->voxels;
	delete[] this->activeTopLevelVoxels;
	delete this->voxelNeighbourGrid;

	int counterSize = 1 << this->dimension;
	for(int i=0; i<counterSize; i++){
		delete[] counter[i];
	}
	delete[] this->counter;
}

void VoxelList::disable(){
	this->disabled = true;
}

Voxel* VoxelList::createVoxel(double* center, double vs, int index){
	double d = 0.5*vs;
	for(unsigned char i=0; i< this->dimension; i++){
		if (center[i] - d < 0.0){
			center[i] = d;
		}else if (center[i] + d > this->size){
			center[i] = this->size - d;
		}
	}
	return new Voxel(this->dimension, center, vs, index);
}


void VoxelList::initVoxels(unsigned char dim){
	int n = (int)(this->size/this->voxelSize + 0.5);
	double *da = new double[dim];
	int *in = new int[dim];

	for(unsigned char i=0; i<dim; i++){
		in[i] = 0;
	}

	int i, index = 0;
	do{
		for(unsigned char i=0; i<dim; i++){
			da[i] = this->voxelSize*in[i]; // da point to the "left bottom" corner of a voxel
		}
		i = position2i(da, dim, this->size, this->voxelSize, n);
		if(index!=i){
			std::cout << "VoxelList::initVoxels: Problem: " << index << " != " << i << std::endl;
		}

		this->voxels[index] = this->createVoxel(da, this->voxelSize*this->dxFactor, index);
		this->activeTopLevelVoxels[index] = true;
		index++;
	}while(increment(in, dim, n-1));
	delete[] da;
	delete[] in;
}


void VoxelList::fillNeighbourGrid(){
	this->voxelNeighbourGrid->clear();
	for(int i=0; i<=this->last; i++){
		this->voxelNeighbourGrid->add(this->voxels[i]);
	}
}

std::unordered_set<Positioned *> * VoxelList::getNeighbours(Voxel *v){
	return this->voxelNeighbourGrid->getNeighbours(v->getPosition());
}


void VoxelList::checkIndexes(){
	for(int i=0; i<=this->last; i++){
		if(this->voxels[i]->index!=i)
			std::cout << "VoxelList::checkIndexes: Error " << i << std::endl;
	}
}


void VoxelList::remove(Voxel *v){
	if (this->disabled)
		return;

	double* vpos = v->getPosition();
	int index = position2i(vpos, this->dimension, this->size, this->initialVoxelSize, (int)(this->size/this->initialVoxelSize + 0.5));
	this->activeTopLevelVoxels[index]=false;

	index = v->index;

	if (index!=last){
		Voxel *vl = this->voxels[last];
		vl->index = index;
		this->voxels[index] = vl;
	}
	this->voxels[last] = NULL;
	this->last--;
	this->voxelNeighbourGrid->remove(v);
	delete v;

//		this.checkIndexes();
	}

bool VoxelList::analyzeVoxel(Voxel *v, Shape *s, BoundaryConditions *bc){

	double* vpos = v->getPosition();
	int index = position2i(vpos, this->dimension, this->size, this->initialVoxelSize, (int)(this->size/this->initialVoxelSize + 0.5));
	if(this->activeTopLevelVoxels[index]==false)
		return true;


	double *da = new double[this->dimension];
	int counterSize = 1 << this->dimension;
	for(int i=0; i<counterSize; i++){
		for(unsigned char j=0; j<this->dimension; j++){
			da[j] = vpos[j] + this->counter[i][j]*this->voxelSize;

			if( !(s->pointInside(bc, da)) ){
				delete[] da;
				return false;
			}
		}
	}

	delete[] da;
	return true;

}

// returns true when the whole voxel is inside an exclusion area of one shape, thus, should be removed
bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid *nl, std::unordered_set<Positioned *> *neighbours, BoundaryConditions *bc){

	double* vpos = v->getPosition();

	// checking if initial voxel containing v is active (do not have a shape inside)
	int index = position2i(vpos, this->dimension, this->size, this->initialVoxelSize, (int)(this->size/this->initialVoxelSize + 0.5));
	if(this->activeTopLevelVoxels[index]==false)
		return true;

	std::unordered_set<Positioned *> vAll;
	std::unordered_set<Positioned *>::iterator itN, itA;

	double *da = new double[this->dimension];

	int counterSize = 1 << this->dimension;
	for(int i=0; i<counterSize; i++){
		for(unsigned char j=0; j<this->dimension; j++){
			da[j] = vpos[j] + this->counter[i][j]*this->voxelSize;
		}
		std::unordered_set<Positioned *> *vNeighbours = (neighbours==NULL) ? nl->getNeighbours(da) : neighbours;
		// at the beginning we add all neighbouring shapes
		if (vAll.size()==0){
			vAll.insert(vNeighbours->begin(), vNeighbours->end());
		}else{
			for(itA = vAll.begin(); itA != vAll.end(); ){
				if( vNeighbours->find(*itA) == vNeighbours->end() ){
					itA = vAll.erase(itA);
				}else{
					itA++;
				}
			}
		}
		for(itA = vAll.begin(); itA != vAll.end(); ){
			if(! ((Shape*)(*itA))->pointInside(bc, da)){
				itA = vAll.erase(itA);
			}else{
				itA++;
			}
		}
		if (vAll.size()==0){
			delete[] da;
			return false;
		}
	}
//	v.analyzed = true;
	delete[] da;
	return true;
}

bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid *nl, BoundaryConditions *bc, int timestamp){
	if (v->lastAnalyzed < timestamp && !this->disabled)
		return this->analyzeVoxel(v, nl, NULL, bc);
	return false;
}

bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid *nl, BoundaryConditions *bc){
	return this->analyzeVoxel(v, nl, NULL, bc);
}

bool VoxelList::analyzeVoxel(Voxel *v, std::unordered_set<Positioned *> *neighbours, BoundaryConditions *bc){
	return this->analyzeVoxel(v, NULL, neighbours, bc);
}

bool VoxelList::splitVoxels(double minDx, int maxVoxels, NeighbourGrid *nl, BoundaryConditions *bc){
	if (this->disabled)
		return false;
	if ((this->voxelSize<2*minDx && pow(2, this->dimension)*this->last > this->beginningVoxelNumber) || pow(2, this->dimension)*this->last > maxVoxels){
//		for(int i=0; i<this.last; i++)
//			this.voxels[i].missCounter = 0;
		return false;
	}

	Voxel** newList = new Voxel*[ ((int)round( pow(2, this->dimension)))*(this->last+1) ];
	this->voxelSize = (this->voxelSize/2.0)*this->dxFactor;
	delete this->distribution;
	this->distribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);

	int index = 0;

	int *in = new int[this->dimension];
	double *da = new double[this->dimension];
	for(int i=0; i<=this->last; i++){
		Voxel *v = this->voxels[i];
		double* vpos = v->getPosition();

		for(unsigned char j=0; j < this->dimension; j++){
			in[j] = 0;
			da[j] = vpos[j];
		}
		do{
			for(unsigned char j=0; j < this->dimension; j++){
				da[j] = vpos[j] + in[j]*this->voxelSize;
			}
			Voxel *vTmp = this->createVoxel(da, this->voxelSize, index);
			if (nl==NULL || bc==NULL || !this->analyzeVoxel(vTmp, nl, NULL, bc)){
				newList[index] = vTmp;
				index++;
			}else{
				delete vTmp;
			}
		}while(increment(in, this->dimension, (unsigned char)1));
		delete this->voxels[i];

		if (i%10000 == 0){ std::cout << "."; std::cout.flush(); }
	}
	delete[] in;
	delete[] da;

	delete[] this->voxels;
	this->voxelNeighbourGrid->clear();

	this->last = index-1;
	this->voxels = newList;
	this->fillNeighbourGrid();
//	this->checkIndexes();
	return true;
}

Voxel * VoxelList::getRandomVoxel(RND *rnd){
	double d = rnd->nextValue();
	return this->voxels[(int)(d*(this->last+1))];
}

double * VoxelList::getRandomPosition(double *result, Voxel *v, RND *rnd){
	double *vpos = v->getPosition();
	for (int i=0; i < this->dimension; i++)
		result[i] = vpos[i] + rnd->nextValue(this->distribution);
	return result;
}

double VoxelList::getVoxelSize(){
	return this->voxelSize;
}

Voxel* VoxelList::get(int i){
	return this->voxels[i];
}

int VoxelList::length(){
	return this->last+1;
}

double VoxelList::getVoxelsSurface(){
	return (this->last+1)*pow(this->voxelSize, this->dimension);
}

std::string VoxelList::toPovray(){
	std::string sRes = "";
	for(int i=0; i<=this->last; i++){
		double *da = this->voxels[i]->getPosition();
		double x1 = da[0], x2 = da[0] + this->voxelSize;
		double y1 = da[1], y2 = da[1] + this->voxelSize;

		sRes += "  polygon {4, < " + std::to_string(x1) + ", " + std::to_string(y1) + ", 1.0>, "
				+ "< " + std::to_string(x1) + ", " + std::to_string(y2) + ", 1.0>, "
				+ "< " + std::to_string(x2) + ", " + std::to_string(y2) + ", 1.0>, "
				+ "< " + std::to_string(x2) + ", " + std::to_string(y1) + ", 1.0> "
				+ " texture { finish { ambient 1 diffuse 0 } pigment { color Black} } }\r\n";
	}
	return sRes;
}
