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

VoxelList::VoxelList(int dim, double s, double d){
	this->dimension = dim;
	this->size = s;
	int n = (int)(s/d)+1;
	this->voxelSize = (s/n);

	n = (int)(this->size/this->voxelSize + 0.5);
	int voxelsLength = (int)(pow(n, dim)+0.5);
	this->voxels = new Voxel*[voxelsLength];
	this->voxelNeighbourGrid = new NeighbourGrid(dim, s, this->voxelSize);

	this->da = new double[dim];
	this->in = new int[dim];

	this->initVoxels(dim);
	this->voxelSize *= this->dxFactor;
	this->beginningVoxelNumber = voxelsLength;
	this->last = voxelsLength-1;
	this->fillNeighbourGrid();


//		this.checkIndexes();
	}

VoxelList::~VoxelList() {
	for(int i=0; i<=this->last; i++){
		delete this->voxels[i];
	}
	delete[] this->voxels;
	delete this->voxelNeighbourGrid;

	delete[] this->da;
	delete[] this->in;
}


Voxel* VoxelList::createVoxel(double* center, double vs, int index){
	return new Voxel(this->dimension, center, vs, index);
}


void VoxelList::initVoxels(int N){
	int n = (int)(this->size/this->voxelSize + 0.5);

	for(int i=0; i<N; i++){
		this->in[i] = 0;
	}

	int index = 0;
	do{
		for(int i=0; i<N; i++){
			this->da[i] = this->voxelSize*(in[i] + 0.5);
		}
		if(index!=position2i(this->da, N, this->size, this->voxelSize, n)){
			std::cout << "Problem" << std::endl;
		}

		this->voxels[index] = this->createVoxel(this->da, this->voxelSize*this->dxFactor, index);
		index++;
	}while(increment(this->in, N, n-1));
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
			std::cout << "Error " << i << std::endl;
	}
}


void VoxelList::remove(Voxel *v){
	int index = v->index;

	if (index!=last){
		Voxel *vl = this->voxels[last];
		vl->index = index;
		this->voxels[index] = vl;
	}
	this->voxelNeighbourGrid->remove(v);
	delete v;
	this->voxels[last] = NULL;
	this->last--;
//		this.checkIndexes();
	}

// returns true when the whole voxel is inside an exclusion area of one shape, thus, should be removed
bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid *nl, std::unordered_set<Positioned *> *neighbours, BoundaryConditions *bc){

	std::unordered_set<Positioned *> vAll;
	std::unordered_set<Positioned *>::iterator itN, itA;

	double* vpos = v->getPosition();

	for(int j=0; j<this->dimension; j++){
		this->in[j] = 0;
		this->da[j] = vpos[j];
	}

	do{
		for(int j=0; j<this->dimension; j++){
			this->da[j] = vpos[j] + (this->in[j]-0.5)*this->voxelSize;
		}
		std::unordered_set<Positioned *> *vNeighbours = (neighbours==NULL) ? nl->getNeighbours(this->da) : neighbours;
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
			if(! ((Shape*)(*itA))->pointInside(bc, this->da)){
				itA = vAll.erase(itA);
			}else{
				itA++;
			}
		}
		if (vAll.size()==0){
			return false;
		}
	}while(increment(this->in, this->dimension, 1));
//	v.analyzed = true;
	return true;
}

bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid *nl, BoundaryConditions *bc, int timestamp){
	if (v->lastAnalyzed < timestamp)
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
	if ((this->voxelSize<2*minDx && pow(2, this->dimension)*this->last > this->beginningVoxelNumber) || pow(2, this->dimension)*this->last > maxVoxels){
//		for(int i=0; i<this.last; i++)
//			this.voxels[i].missCounter = 0;
		return false;
	}

	Voxel** newList = new Voxel*[ ((int)round( pow(2, this->dimension)))*(this->last+1) ];
	double newDx = (this->voxelSize/2.0)*this->dxFactor;
	int index = 0;

	double *da = new double[this->dimension];
	int *in = new int[this->dimension];

	for(int i=0; i<=this->last; i++){
		Voxel *v = this->voxels[i];
		double* vpos = v->getPosition();


		for(int j=0; j < this->dimension; j++){
			in[j] = 0;
			da[j] = vpos[j];
		}
		do{
			for(int j=0; j < this->dimension; j++){
				da[j] = vpos[j] + (in[j]-0.5)*newDx;
			}
			Voxel *vTmp = new Voxel(this->dimension, da, newDx, index);
			if (nl==NULL || bc==NULL || !this->analyzeVoxel(vTmp, nl, NULL, bc)){
				newList[index] = vTmp;
				index++;
			}else{
				delete vTmp;
			}
		}while(increment(in, this->dimension, 1));
		delete this->voxels[i];
	}

	delete[] da;
	delete[] in;

	delete[] this->voxels;
	this->voxelNeighbourGrid->clear();

	this->last = index-1;
	this->voxels = newList;
	this->voxelSize = newDx;
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
		result[i] = vpos[i] + (rnd->nextValue()-0.5)*this->voxelSize;
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

/*
	synchronized public void drawVoxels(Graphics g, double scale, double[] ta){
		int w = (int)(this.voxelSize*scale);
		w += 1;
		for(int i=0; i<=this.last; i++){
			Voxel v = this.voxels[i];
			g.fillRect((int)((v.center[0]-0.5*this.voxelSize+ta[0])*scale), (int)((v.center[1]-0.5*this.voxelSize+ta[1])*scale), w, w);
		}
	}
}
*/
