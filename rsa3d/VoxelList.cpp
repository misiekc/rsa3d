/*
 * VoxelList.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */


#include "VoxelList.h"

#include <iostream>
#include <cmath>

#include "Voxel.h"
#include "math.h"
#include "NeighbourGrid.h"

#include "Utils.cpp"

VoxelList::VoxelList(int dim, double s, double d){
	this->dimension = dim;
	this->size = s;
	int n = (int)(s/d)+1;
	this->voxelSize = (s/n);

	n = (int)(this->size/this->voxelSize + 0.5);
	int voxelsLength = (int)(pow(n, dim)+0.5);
	this->voxels = new Voxel[voxelsLength];
	this->voxelNeighbourGrid = new NeighbourGrid(dim, s, this->voxelSize);

	this->initVoxels(dim);
	this->voxelSize *= dxFactor;
	this->beginningVoxelNumber = voxelsLength;
	this->last = voxelsLength-1;
	this->fillNeighbourGrid();

//		this.checkIndexes();
	}

VoxelList::~VoxelList() {
	// TODO Auto-generated destructor stub
}


Voxel* VoxelList::createVoxel(double* center, double vs, int index){
	return new Voxel(this->dimension, center, vs, index);
}


void VoxelList::initVoxels(int N){
	int n = (int)(this->size/this->voxelSize + 0.5);
	double da[N];
	int in[N];

	for(int i=0; i<N; i++){
		in[i] = 0;
	}

	int index = 0;
	do{
		for(int i=0; i<N; i++){
			da[i] = this->voxelSize*(in[i] + 0.5);
		}
		if(index!=position2i(da, N, this->size, this->voxelSize, n)){
			std::cout << "Problem" << std::endl;
		}

		this->voxels[index] = this->createVoxel(da, this->voxelSize*dxFactor, index);
		index++;
	}while(increment(in, N, n-1));
}


void VoxelList::fillNeighbourGrid(){
	this->voxelNeighbourGrid->clear();
	for(int i=0; i<=this->last; i++){
		this->voxelNeighbourGrid->add(this->voxels + i);
	}
}

Voxel* VoxelList::getNeighbours(Voxel *v){
	return this->voxelNeighbourGrid->getNeighbours(v->getPosition());
}


void VoxelList::checkIndexes(){
	for(int i=0; i<=this->last; i++){
		if(this->voxels[i].index!=i)
			std::cout << "Error " << i << std::endl;
	}
}


void VoxelList::remove(Voxel *v){
	int index = v->index;

	if (index!=last){
		Voxel *vl = this->voxels + last;
		vl->index = index;
		this->voxels[index] = vl;
	}
	this->voxels[last] = NULL;
	this->last--;
	this->voxelNeighbourGrid->remove(v);
	delete v;

//		this.checkIndexes();
	}

bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid *nl, std::vector<Shape*> *neighbours, BoundaryConditions *bc){

	std::vector<Shape *> vAll;
	std::vector<Shape *>::iterator it1, it2;

	int in[this->dimension];
	double da[this->dimension];
	double* vpos = v->getPosition();

	for(int j=0; j<this->dimension; j++){
		in[j] = 0;
		da[j] = vpos[j];
	}

	do{
		for(int j=0; j<this->dimension; j++){
			da[j] = vpos[j] + (in[j]-0.5)*this->voxelSize;
		}
		std::vector<Shape*> *vNeighbours = (neighbours==NULL) ? nl->getNeighbours(da) : neighbours;
		if (vAll.size()==0){
			vAll.insert(vAll.end(), vNeighbours->begin(), vNeighbours->end());
		}else{
			for(it1 = vAll.begin(); it1 != vAll.end(); it1++){
				if( (it2 = std::find(vNeighbours->begin(), vNeighbours->end(), *it1)) != vNeighbours->end()){
					vAll.erase(it1);
					it1--;
				}
			}
		}
		for(it1 = vAll.begin(); it1 != vAll.end(); it1++){
			if(! (*it1)->pointInside(bc, da)){
				vAll.erase(it1);
				it1--;
			}
		}
	}while(increment(in, this->dimension, 1));
//	v.analyzed = true;
	if (vAll.size()>0){
		return true;
	}else{
		return false;
	}
}

bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid *nl, BoundaryConditions *bc, int timestamp){
	if (v->lastAnalyzed < timestamp)
		return this->analyzeVoxel(v, nl, NULL, bc);
	return false;
}

bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid *nl, BoundaryConditions *bc){
	return this->analyzeVoxel(v, nl, NULL, bc);
}

bool VoxelList::analyzeVoxel(Voxel *v, std::vector<Shape *> *neighbours, BoundaryConditions *bc){
	return this->analyzeVoxel(v, NULL, neighbours, bc);
}

/*
	public boolean splitVoxels(double minDx, int maxVoxels){
		if ((this.voxelSize<2*minDx && Math.pow(2, this.dimension)*this.last>this.beginningVoxelsNumber) || Math.pow(2, this.dimension)*this.last>maxVoxels){
//			for(int i=0; i<this.last; i++)
//				this.voxels[i].missCounter = 0;
			return false;
		}
		Commons.logger.info("Voxel size " + this.voxelSize);
		Voxel[] newList = new Voxel[ ((int)Math.round(Math.pow(2, this.dimension)))*(this.last+1) ];
		double newDx = (this.voxelSize/2.0)*dxFactor;
		int index = 0;
		for(int i=0; i<=this.last; i++){
			Voxel v = this.voxels[i];

			int[] in = new int[v.center.length];
			double[] da = v.center.clone();
			for(int j=0; j<in.length; j++){
				in[j] = 0;
			}
			do{
				for(int j=0; j<da.length; j++){
					da[j] = v.center[j] + (in[j]-0.5)*newDx;
				}
				newList[index] = new Voxel(da, newDx, index);
				index++;
			}while(Commons.increment(in, 1));
			this.voxels[i] = null;
		}
		this.last = newList.length-1;
		this.voxels = newList;
		this.voxelSize = newDx;
//		this.checkIndexes();
		Commons.logger.info("Voxel size " + this.voxelSize);
		return true;
	}
*/

bool VoxelList::splitVoxels(double minDx, int maxVoxels, NeighbourGrid *nl, BoundaryConditions *bc){
	if ((this->voxelSize<2*minDx && pow(2, this->dimension)*this->last > this->beginningVoxelNumber) || pow(2, this->dimension)*this->last > maxVoxels){
//		for(int i=0; i<this.last; i++)
//			this.voxels[i].missCounter = 0;
		return false;
	}

	Voxel* newList = new Voxel[ ((int)round( pow(2, this->dimension)))*(this->last+1) ];
	double newDx = (this->voxelSize/2.0)*VoxelList::dxFactor;
	int index = 0;
	for(int i=0; i<=this->last; i++){
		Voxel v = this->voxels[i];
		double* vpos = v.getPosition();

		int in[this->dimension];
		double da[this->dimension];
		for(int j=0; j < this->dimension; j++){
			in[j] = 0;
			da[i] = vpos[i];
		}
		do{
			for(int j=0; j < this->dimension; j++){
				da[j] = vpos[j] + (in[j]-0.5)*newDx;
			}
			Voxel *vTmp = new Voxel(this->dimension, da, newDx, index);
			if (nl==NULL || bc==NULL || !this->analyzeVoxel(vTmp, nl, NULL, bc)){
				newList[index] = vTmp;
				index++;
			}
		}while(increment(in, this->dimension, 1));
		delete this->voxels[i];
	}
	this->last = index-1;
	this->voxels = newList;
	this->voxelSize = newDx;
	this->fillNeighbourList();
//	this->checkIndexes();
	return true;
}

Voxel * VoxelList::getRandomVoxel(RND *rnd){
	return &(this->voxels[(int)(rnd->nextValue()*(this->last+1)]);
}

double * VoxelList::getRandomPosition(double *result, Voxel *v, RND *rnd){
	double da[this->dimension];
	double *vpos = v->getPosition();
	for (int i=0; i < this->dimension; i++)
		da[i] = vpos[i] + (rnd->nextValue()-0.5)*this->voxelSize;
	return da;
}

double VoxelList::getVoxelSize(){
	return this->voxelSize;
}

Voxel* VoxelList::get(int i){
	return this->voxels + i;
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
