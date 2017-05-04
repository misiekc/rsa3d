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

VoxelList::VoxelList(unsigned char dim, double s, double d){
	this->dimension = dim;
	this->size = s;
	int n = (int)(s/d)+1;
	this->voxelSize = (s/n);

	n = (int)(this->size/this->voxelSize + 0.5);
	int voxelsLength = (int)(pow(n, dim)+0.5);
	this->voxels = new Voxel*[voxelsLength];
	this->voxelNeighbourGrid = new NeighbourGrid(dim, s, this->voxelSize);

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
}


Voxel* VoxelList::createVoxel(double* center, double vs, int index){
	return new Voxel(this->dimension, center, vs, index);
}


void VoxelList::initVoxels(unsigned char dim){
	int n = (int)(this->size/this->voxelSize + 0.5);
	double *da = new double[dim];
	int *in = new int[dim];

	for(unsigned char i=0; i<dim; i++){
		in[i] = 0;
	}

	int index = 0;
	do{
		for(unsigned char i=0; i<dim; i++){
			da[i] = this->voxelSize*(in[i] + 0.5);
		}
		if(index!=position2i(da, dim, this->size, this->voxelSize, n)){
			std::cout << "Problem" << std::endl;
		}

		this->voxels[index] = this->createVoxel(da, this->voxelSize*this->dxFactor, index);
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
	this->voxels[last] = NULL;
	this->last--;
	this->voxelNeighbourGrid->remove(v);
	delete v;

//		this.checkIndexes();
	}

// returns true when the whole voxel is inside an exclusion area of one shape, thus, should be removed
bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid *nl, std::unordered_set<Positioned *> *neighbours, BoundaryConditions *bc){

	std::unordered_set<Positioned *> vAll;
	std::unordered_set<Positioned *>::iterator itN, itA;

	int *in = new int[this->dimension];
	double *da = new double[this->dimension];
	double* vpos = v->getPosition();

	for(unsigned char j=0; j<this->dimension; j++){
		in[j] = 0;
		da[j] = vpos[j];
	}

	do{
		for(unsigned char j=0; j<this->dimension; j++){
			da[j] = vpos[j] + (in[j]-0.5)*this->voxelSize;
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
			delete[] in;
			delete[] da;
			return false;
		}
	}while(increment(in, this->dimension, 1));
//	v.analyzed = true;
	delete[] in;
	delete[] da;
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
				da[j] = vpos[j] + (in[j]-0.5)*newDx;
			}
			Voxel *vTmp = new Voxel(this->dimension, da, newDx, index);
			if (nl==NULL || bc==NULL || !this->analyzeVoxel(vTmp, nl, NULL, bc)){
				newList[index] = vTmp;
				index++;
			}else{
				delete vTmp;
			}
		}while(increment(in, this->dimension, (unsigned char)1));
		delete this->voxels[i];
	}
	delete[] in;
	delete[] da;

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
