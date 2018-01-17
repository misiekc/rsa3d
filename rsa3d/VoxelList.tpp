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

/**
 * d - requested initial size of a voxel
 */
template <ushort DIMENSION>
double VoxelList<DIMENSION>::findInitialVoxelSize(double d){
	double dRes = 1.0;
	while (dRes > d)
		dRes /= 2;
	while (2*dRes < d)
		dRes *= 2;
	return dRes;
}

template <ushort DIMENSION>
inline int VoxelList<DIMENSION>::getLinearNumberOfVoxels(double vs){
	return (int)(this->size/vs) + 1;
}

/**
dim - packing dimension
s - packing size (linear)
d - requested initial size of a voxel
**/
template <ushort DIMENSION>
VoxelList<DIMENSION>::VoxelList(double s, double d){
	this->size = s;
	this->voxelSize = this->findInitialVoxelSize(d);
	this->initialVoxelSize = this->voxelSize;

	this->distribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);
	this->disabled = false;

	int n = this->getLinearNumberOfVoxels(this->voxelSize);
	int voxelsLength = (int)(pow(n, DIMENSION)+0.5);
	this->voxels = new Voxel<DIMENSION>*[voxelsLength];
	this->activeTopLevelVoxels = new bool[voxelsLength];
	this->voxelNeighbourGrid = new NeighbourGrid<Voxel<DIMENSION>>(DIMENSION, this->voxelSize*n, n);


	int in[DIMENSION];
	for(ushort i=0; i<DIMENSION; i++)
		in[i] = 0;
	int index = 0;
	do{
		std::copy(in, in+DIMENSION, offset[index]);
		index++;
	}while(increment(in, DIMENSION, 1));

	this->initVoxels();
	this->voxelSize *= this->dxFactor;
	this->beginningVoxelNumber = voxelsLength;
	this->last = voxelsLength-1;
	this->fillNeighbourGrid();

//		this.checkIndexes();
	}

template <ushort DIMENSION>
VoxelList<DIMENSION>::~VoxelList() {
	delete this->distribution;

	for(int i=0; i<=this->last; i++){
		delete this->voxels[i];
	}
	delete[] this->voxels;
	delete[] this->activeTopLevelVoxels;
	delete this->voxelNeighbourGrid;
}

template <ushort DIMENSION>
void VoxelList<DIMENSION>::disable(){
	this->disabled = true;
}

template <ushort DIMENSION>
Voxel<DIMENSION>* VoxelList<DIMENSION>::createVoxel(double* leftbottom, double vs, int index){
/*
	for(int i=0; i< this->dimension; i++){
		if (leftbottom[i] < 0.0){
			leftbottom[i] = 0.0;
		}else if (leftbottom[i] + vs > this->size){
			leftbottom[i] = this->size - vs;
		}
	}
*/
	return new Voxel<DIMENSION>(leftbottom, vs, index);
}

template <ushort DIMENSION>
void VoxelList<DIMENSION>::initVoxels(){
	int n = this->getLinearNumberOfVoxels(this->voxelSize);
	double da[DIMENSION];
	int in[DIMENSION];

	for(ushort i = 0; i<DIMENSION; i++)
		in[i] = 0;

	int i, index = 0;
	do{
		for(unsigned char i=0; i<DIMENSION; i++){
			da[i] = this->voxelSize*in[i]; // da point to the "left bottom" corner of a voxel
		}
		i = position2i(da, DIMENSION, n*this->voxelSize, this->voxelSize, n);
		if(index!=i){
			std::cout << "VoxelList::initVoxels: Problem: " << index << " != " << i << std::endl;
		}

		this->voxels[index] = this->createVoxel(da, this->voxelSize*this->dxFactor, index);
		this->activeTopLevelVoxels[index] = true;
		index++;
	}while(increment(in, DIMENSION, n-1));
}

template <ushort DIMENSION>
void VoxelList<DIMENSION>::fillNeighbourGrid(){
	this->voxelNeighbourGrid->clear();
	for(int i=0; i<=this->last; i++){
		this->voxelNeighbourGrid->add(this->voxels[i], this->voxels[i]->getPosition());
	}
}

template <ushort DIMENSION>
void VoxelList<DIMENSION>::getNeighbours(std::vector<Voxel<DIMENSION> *> *result, Voxel<DIMENSION> *v){
	return this->voxelNeighbourGrid->getNeighbours(result, v->getPosition());
}


template <ushort DIMENSION>
void VoxelList<DIMENSION>::checkIndexes(){
	for(int i=0; i<=this->last; i++){
		if(this->voxels[i]->index!=i)
			std::cout << "VoxelList::checkIndexes: Error " << i << std::endl;
	}
}

template <ushort DIMENSION>
void VoxelList<DIMENSION>::remove(Voxel<DIMENSION> *v){
	if (this->disabled)
		return;

	int index = v->index;

	if (index!=last){
		Voxel<DIMENSION> *vl = this->voxels[last];
		vl->index = index;
		this->voxels[index] = vl;
	}
	this->voxels[last] = NULL;
	this->last--;
	this->voxelNeighbourGrid->remove(v, v->getPosition());
	delete v;

//		this.checkIndexes();
}

template <ushort DIMENSION>
void VoxelList<DIMENSION>::removeTopLevelVoxel(Voxel<DIMENSION> *v){
	if (this->disabled)
		return;

	double* vpos = v->getPosition();
	int n = (int)(this->size/this->initialVoxelSize) + 1;
	int index = position2i(vpos, DIMENSION, n*this->initialVoxelSize, this->initialVoxelSize, n);
	this->activeTopLevelVoxels[index]=false;
}


template <ushort DIMENSION>
Voxel<DIMENSION> * VoxelList<DIMENSION>::getVoxel(double* da){
	std::vector<Voxel<DIMENSION> *> *vTmp = this->voxelNeighbourGrid->getCell(da);
	for(Voxel<DIMENSION> *v : *vTmp){
		if (v->isInside(da, this->voxelSize)){
			return v;
		}
	}
	return NULL;
}

template <ushort DIMENSION>
bool VoxelList<DIMENSION>::analyzeVoxel(Voxel<DIMENSION> *v, Shape<DIMENSION> *s, BoundaryConditions *bc){

	double* vpos = v->getPosition();
	int n = (int)(this->size/this->initialVoxelSize) + 1;
	int index = position2i(vpos, DIMENSION, n*this->initialVoxelSize, this->initialVoxelSize, n);
	if(this->activeTopLevelVoxels[index]==false)
		return true;


	double da[DIMENSION];
	int counterSize = 1 << DIMENSION;
	for(int i=0; i<counterSize; i++){
		for(ushort j=0; j<DIMENSION; j++){
			da[j] = vpos[j] + this->offset[i][j]*this->voxelSize;

			if( !(s->pointInside(bc, da)) ){
				return false;
			}
		}
	}
	return true;
}

/*
template <ushort DIMENSION>
bool VoxelList<DIMENSION>::analyzeVoxel(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, std::unordered_set<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc){
	bool b1 = analyzeVoxelOLD(v, nl, neighbours, bc);
	bool b2 = analyzeVoxelNEW(v, nl, neighbours, bc);
	if (b1!=b2){
		double *vPos = v->getPosition();
		std::cout << "PROBLEM: (" << vPos[0] << ", " << vPos[1] << ", " << vPos[2] << ")" << std::endl;
	}
	return b1;
}
// returns true when the whole voxel is inside an exclusion area of one shape, thus, should be removed
template <ushort DIMENSION>
bool VoxelList<DIMENSION>::analyzeVoxelOLD(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, std::unordered_set<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc){

	double* vpos = v->getPosition();

	// checking if initial voxel containing v is active (do not have a shape inside)
	int n = this->getLinearNumberOfVoxels(this->initialVoxelSize);
	int index = position2i(vpos, DIMENSION, n*this->initialVoxelSize, this->initialVoxelSize, n);
	if(this->activeTopLevelVoxels[index]==false)
		return true;

	std::unordered_set<Shape<DIMENSION> *> vAll;
	typename std::unordered_set<Shape<DIMENSION> *>::iterator itN, itA;
	std::unordered_set<Shape<DIMENSION> *> *vNeighbours;

	if(neighbours==NULL){
		vNeighbours = new std::unordered_set<Shape<DIMENSION> *>();
	}else{
		vNeighbours = neighbours;
	}

	double da[DIMENSION];

	int counterSize = 1 << DIMENSION;
	for(int i=0; i<counterSize; i++){
		for(ushort j=0; j<DIMENSION; j++){
			da[j] = vpos[j] + this->offset[i][j]*this->voxelSize;
		}

		if (neighbours==NULL){
			nl->getNeighbours(vNeighbours, da);
		}
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
			if(! (*itA)->pointInside(bc, da)){
				itA = vAll.erase(itA);
			}else{
				itA++;
			}
		}
		if (vAll.size()==0){
			if (neighbours==NULL)
				delete vNeighbours;
			return false;
		}
	}
//	v.analyzed = true;
	if (neighbours==NULL)
		delete vNeighbours;
	return true;
}
*/

// returns true when the whole voxel is inside an exclusion area of one shape, thus, should be removed
template <ushort DIMENSION>
bool VoxelList<DIMENSION>::analyzeVoxel(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, std::vector<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc){

	double* vpos = v->getPosition();

	// checking if initial voxel containing v is active (do not have a shape inside)
	int n = this->getLinearNumberOfVoxels(this->initialVoxelSize);
	int index = position2i(vpos, DIMENSION, n*this->initialVoxelSize, this->initialVoxelSize, n);
	if(this->activeTopLevelVoxels[index]==false)
		return true;

//	std::unordered_set<Shape<DIMENSION> *> vAll;
//	typename std::unordered_set<Shape<DIMENSION> *>::iterator itN, itA;
	std::vector<Shape<DIMENSION> *> *vNeighbours;

	if(neighbours==NULL){
		vNeighbours = new std::vector<Shape<DIMENSION> *>();
		nl->getNeighbours(vNeighbours, vpos);
	}else{
		vNeighbours = neighbours;
	}

	double da[DIMENSION];
	int counterSize = 1 << DIMENSION;
	bool shapeCoversVertices = false;

	for(Shape<DIMENSION> *s : *vNeighbours){
		shapeCoversVertices = true;
		for(int i=0; i<counterSize; i++){
			for(ushort j=0; j<DIMENSION; j++){
				da[j] = vpos[j] + this->offset[i][j]*this->voxelSize;
			}
			if(!s->pointInside(bc, da)){
				shapeCoversVertices = false;
				break;
			}
		}
		if (shapeCoversVertices)
			break;
	}
	if (neighbours==NULL)
		delete vNeighbours;
	return shapeCoversVertices;
}


template <ushort DIMENSION>
bool VoxelList<DIMENSION>::analyzeVoxel(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, BoundaryConditions *bc, int timestamp){
	if (v->lastAnalyzed < timestamp && !this->disabled)
		return this->analyzeVoxel(v, nl, NULL, bc);
	return false;
}

template <ushort DIMENSION>
bool VoxelList<DIMENSION>::analyzeVoxel(Voxel<DIMENSION> *v, NeighbourGrid<Shape<DIMENSION>> *nl, BoundaryConditions *bc){
	return this->analyzeVoxel(v, nl, NULL, bc);
}

template <ushort DIMENSION>
bool VoxelList<DIMENSION>::analyzeVoxel(Voxel<DIMENSION> *v, std::vector<Shape<DIMENSION> *> *neighbours, BoundaryConditions *bc){
	return this->analyzeVoxel(v, NULL, neighbours, bc);
}

template <ushort DIMENSION>
bool VoxelList<DIMENSION>::splitVoxels(double minDx, int maxVoxels, NeighbourGrid<Shape<DIMENSION>> *nl, BoundaryConditions *bc){
	if (this->disabled)
		return false;
	if ((this->voxelSize<2*minDx && pow(2, DIMENSION)*this->last > this->beginningVoxelNumber) || pow(2, DIMENSION)*this->last > maxVoxels){
//		for(int i=0; i<this.last; i++)
//			this.voxels[i].missCounter = 0;
		return false;
	}

	Voxel<DIMENSION>** newList = new Voxel<DIMENSION>*[ ((int)round( pow(2, DIMENSION)))*(this->last+1) ];
	this->voxelSize = (this->voxelSize/2.0)*this->dxFactor;
	delete this->distribution;
	this->distribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);

	int index = 0;

	int in[DIMENSION];
	double da[DIMENSION];
	for(int i=0; i<=this->last; i++){
		Voxel<DIMENSION> *v = this->voxels[i];
		double* vpos = v->getPosition();

		for(ushort j=0; j < DIMENSION; j++){
			in[j] = 0;
			da[j] = vpos[j];
		}
		do{
			bool doCreate = true;
			for(ushort j=0; j < DIMENSION; j++){
				da[j] = vpos[j] + in[j]*this->voxelSize;
				if (da[j]>this->size)
					doCreate = false;
			}
			if (doCreate){
				Voxel<DIMENSION> *vTmp = this->createVoxel(da, this->voxelSize, index);
				if (nl==NULL || bc==NULL || !this->analyzeVoxel(vTmp, nl, NULL, bc)){
					newList[index] = vTmp;
					index++;
				}else{
					delete vTmp;
				}
			}
		}while(increment(in, DIMENSION, (unsigned char)1));
		delete this->voxels[i];

		if (i%10000 == 0){ std::cout << "."; std::cout.flush(); }
	}

	delete[] this->voxels;
	this->voxelNeighbourGrid->clear();

	this->last = index-1;
	this->voxels = newList;
	this->fillNeighbourGrid();
//	this->checkIndexes();
	return true;
}

template <ushort DIMENSION>
Voxel<DIMENSION> * VoxelList<DIMENSION>::getRandomVoxel(RND *rnd){
	double d = rnd->nextValue();
	return this->voxels[(int)(d*(this->last+1))];
}

template <ushort DIMENSION>
double * VoxelList<DIMENSION>::getRandomPosition(double *result, Voxel<DIMENSION> *v, RND *rnd){
	double *vpos = v->getPosition();
	for (ushort i=0; i < DIMENSION; i++)
		result[i] = vpos[i] + rnd->nextValue(this->distribution);
	return result;
}

template <ushort DIMENSION>
double VoxelList<DIMENSION>::getVoxelSize(){
	return this->voxelSize;
}

template <ushort DIMENSION>
Voxel<DIMENSION>* VoxelList<DIMENSION>::get(int i){
	return this->voxels[i];
}

template <ushort DIMENSION>
int VoxelList<DIMENSION>::length(){
	return this->last+1;
}

template <ushort DIMENSION>
double VoxelList<DIMENSION>::getVoxelsSurface(){
	return (this->last+1)*pow(this->voxelSize, DIMENSION);
}

template <ushort DIMENSION>
std::string VoxelList<DIMENSION>::toPovray(){
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
/*
		sRes += "  text { ttf \"timrom.ttf\" \"" + std::to_string(i) + "\" 1, 0 pigment { color White } scale 0.1 translate < ";
		for(unsigned char j=0; j<this->dimension; j++)
			sRes += std::to_string(da[j]+0.5*this->voxelSize) + ", ";
		sRes +=  "1.0003> }\r\n";
*/
	}
	return sRes;
}
