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
inline int VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getLinearNumberOfVoxels(double vs){
	return (int)(this->size/vs) + 1;
}

/**
dim - packing dimension
s - packing size (linear)
d - requested initial size of a voxel
**/
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::VoxelList(double s, double d){
	this->size = s;
	this->voxelSize = this->findInitialVoxelSize(d);
	this->initialVoxelSize = this->voxelSize;
	this->angularSize = 2*M_PI;

	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularSize);
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
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>* VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::createVoxel(double* leftbottom, double* orientation, int index){
/*
	for(int i=0; i< this->dimension; i++){
		if (leftbottom[i] < 0.0){
			leftbottom[i] = 0.0;
		}else if (leftbottom[i] + vs > this->size){
			leftbottom[i] = this->size - vs;
		}
	}
*/
	return new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>(leftbottom, orientation, index);
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

		this->voxels[index] = this->createVoxel(position, orientation, index);
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
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getNeighbours(std::vector<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *result, Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v){
	return this->voxelNeighbourGrid->getNeighbours(result, v->getPosition());
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
void VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::removeTopLevelVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v){
	if (this->disabled)
		return;

	double* vpos = v->getPosition();
	int n = (int)(this->size/this->initialVoxelSize) + 1;
	int index = position2i(vpos, SPATIAL_DIMENSION, n*this->initialVoxelSize, this->initialVoxelSize, n);
	this->activeTopLevelVoxels[index]=false;
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> * VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVoxel(double *pos, double *angle){
	std::vector<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *vTmp = this->voxelNeighbourGrid->getCell(pos);
	for(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v : *vTmp){
		if (v->isInside(pos, this->voxelSize, angle, this->angularSize)){
			return v;
		}
	}
	return NULL;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s, BoundaryConditions *bc){

	double* vpos = v->getPosition();
	int n = (int)(this->size/this->initialVoxelSize) + 1;
	int index = position2i(vpos, SPATIAL_DIMENSION, n*this->initialVoxelSize, this->initialVoxelSize, n);
	if(this->activeTopLevelVoxels[index]==false)
		return true;


	double position[SPATIAL_DIMENSION];
	int counterSize = 1 << SPATIAL_DIMENSION;
	for(int i=0; i<counterSize; i++){
		for(ushort j=0; j<SPATIAL_DIMENSION; j++){
			position[j] = vpos[j] + this->offset[i][j]*this->voxelSize;
			if( !(s->pointInside(bc, position, v->getOrientation(), this->angularSize)) ){
				return false;
			}
		}
	}
	return true;
}

// returns true when the whole voxel is inside an exclusion area of one shape, thus, should be removed
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *neighbours, BoundaryConditions *bc){

	double* vpos = v->getPosition();

	// checking if initial voxel containing v is active (do not have a shape inside)
	int n = this->getLinearNumberOfVoxels(this->initialVoxelSize);
	int index = position2i(vpos, SPATIAL_DIMENSION, n*this->initialVoxelSize, this->initialVoxelSize, n);
	if(this->activeTopLevelVoxels[index]==false)
		return true;

	std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *vNeighbours;

	if(neighbours==NULL){
		vNeighbours = new std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *>();
		nl->getNeighbours(vNeighbours, vpos);
	}else{
		vNeighbours = neighbours;
	}

	double position[SPATIAL_DIMENSION];
	int counterSize = 1 << SPATIAL_DIMENSION;
	bool shapeCoversVertices = false;

	for(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s : *vNeighbours){
		shapeCoversVertices = true;
		for(int i=0; i<counterSize; i++){
			for(ushort j=0; j<SPATIAL_DIMENSION; j++){
				position[j] = vpos[j] + this->offset[i][j]*this->voxelSize;
			}
			if(!s->pointInside(bc, position, v->getOrientation(), this->angularSize)){
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

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc, int timestamp){
	if (v->lastAnalyzed < timestamp && !this->disabled)
		return this->analyzeVoxel(v, nl, NULL, bc);
	return false;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc){
	return this->analyzeVoxel(v, nl, NULL, bc);
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxel(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v, std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *neighbours, BoundaryConditions *bc){
	return this->analyzeVoxel(v, NULL, neighbours, bc);
}

#ifdef _OPENMP

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::splitVoxels(double minDx, int maxVoxels, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc){
	if (this->disabled)
		return false;
	if ((this->voxelSize<2*minDx && pow(2, SPATIAL_DIMENSION+ANGULAR_DIMENSION)*this->last > this->beginningVoxelNumber) || pow(2, SPATIAL_DIMENSION+ANGULAR_DIMENSION)*this->last > maxVoxels){
//		for(int i=0; i<this.last; i++)
//			this.voxels[i].missCounter = 0;
		return false;
	}

	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>** newList = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>*[ ((int)round( pow(2, SPATIAL_DIMENSION+ANGULAR_DIMENSION)))*(this->last+1) ];
	this->voxelSize = (this->voxelSize/2.0)*this->dxFactor;
	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);

	this->angularSize = this->angularSize/2.0;
	delete this->angularDistribution;
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularSize);


	unsigned short spatialLoop = 1 << SPATIAL_DIMENSION;
	unsigned short angularLoop = 1 << ANGULAR_DIMENSION;

	#pragma omp parallel for
	for(int i=0; i<=this->last; i++){

		int inpos[SPATIAL_DIMENSION];
		double position[SPATIAL_DIMENSION];
		int inangle[ANGULAR_DIMENSION];
		double orientation[ANGULAR_DIMENSION];

		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v = this->voxels[i];

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

		for(int j=0; j<spatialLoop; j++){
			bool doCreate = true;
			// checking if a new voxel will be inside a packing
			for(ushort j=0; j < SPATIAL_DIMENSION; j++){
				position[j] = vpos[j] + inpos[j]*this->voxelSize;
				if (position[j]>this->size)
					doCreate = false;
			}
			if (doCreate){
				int k = 0;
				for(ushort j=0; j < ANGULAR_DIMENSION; j++){
					inangle[j] = 0;
				}
				do{
					for(unsigned short j=0; j<ANGULAR_DIMENSION; j++){
						orientation[j] = vangle[j] + inangle[j]*this->angularSize;
					}

					Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *vTmp = this->createVoxel(position, orientation, i*spatialLoop*angularLoop + j*angularLoop + k);
					if (nl==NULL || bc==NULL || !this->analyzeVoxel(vTmp, nl, NULL, bc)){
						newList[i*spatialLoop*angularLoop + j*angularLoop + k] = vTmp;
					}else{
						delete vTmp;
						newList[i*spatialLoop*angularLoop + j*angularLoop + k] = NULL;
					}
					k++;
				}while(increment(inangle, ANGULAR_DIMENSION, (unsigned char)1));
			}else{
				for (int k=0; k<angularLoop; k++)
					newList[i*spatialLoop*angularLoop + j*angularLoop + k] = NULL;
			}
			increment(inpos, SPATIAL_DIMENSION, (unsigned char)1);
		}
		delete this->voxels[i];

		if (i%10000 == 0){ std::cout << "."; std::cout.flush(); }
	}

	delete[] this->voxels;

	int endIndex = this->last*spatialLoop*angularLoop - 1;
	int beginIndex = 0;

	while(newList[endIndex]==NULL)
		endIndex--;

	while (beginIndex<endIndex){
		if(newList[beginIndex]==NULL){
			newList[beginIndex] = newList[endIndex];
			newList[beginIndex]->index = beginIndex;
			newList[endIndex] = NULL;
		}
		beginIndex++;
		while(newList[endIndex]==NULL)
			endIndex--;
	}

	this->voxelNeighbourGrid->clear();

	this->last = endIndex;
	this->voxels = newList;
	this->fillNeighbourGrid();
	this->checkIndexes();
	return true;
}

#else

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::splitVoxels(double minDx, int maxVoxels, NeighbourGrid<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>> *nl, BoundaryConditions *bc){
	if (this->disabled)
		return false;
	if ((this->voxelSize<2*minDx && pow(2, SPATIAL_DIMENSION+ANGULAR_DIMENSION)*this->last > this->beginningVoxelNumber) || pow(2, SPATIAL_DIMENSION+ANGULAR_DIMENSION)*this->last > maxVoxels){
//		for(int i=0; i<this.last; i++)
//			this.voxels[i].missCounter = 0;
		return false;
	}

	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>** newList = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>*[ ((int)round( pow(2, SPATIAL_DIMENSION+ANGULAR_DIMENSION)))*(this->last+1) ];
	this->voxelSize = (this->voxelSize/2.0)*this->dxFactor;
	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->voxelSize);

	this->angularSize = this->angularSize/2.0;
	delete this->angularDistribution;
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularSize);

	int index = 0;

	int inpos[SPATIAL_DIMENSION];
	double position[SPATIAL_DIMENSION];
	int inangle[ANGULAR_DIMENSION];
	double orientation[ANGULAR_DIMENSION];

	for(int i=0; i<=this->last; i++){
		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v = this->voxels[i];

		double* vpos = v->getPosition();
		for(ushort j=0; j < SPATIAL_DIMENSION; j++){
			inpos[j] = 0;
			position[j] = vpos[j];
		}

		double* vangle = v->getOrientation();
		for(ushort j=0; j < ANGULAR_DIMENSION; j++){
			orientation[j] = vangle[j];
		}

		do{
			bool doCreate = true;
			// checking if a new voxel will be inside a packing. it is possible because at the beginning voxels covers more than a packing (due to initial size optimisation)
			for(ushort j=0; j < SPATIAL_DIMENSION; j++){
				position[j] = vpos[j] + inpos[j]*this->voxelSize;
				if (position[j]>this->size)
					doCreate = false;
			}
			if (doCreate){
				for(ushort j=0; j < ANGULAR_DIMENSION; j++){
					inangle[j] = 0;
				}
				do{
					for(unsigned short j=0; j<ANGULAR_DIMENSION; j++){
						orientation[j] = vangle[j] + inangle[j]*this->angularSize;
					}

					Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *vTmp = this->createVoxel(position, orientation, index);
					if (nl==NULL || bc==NULL || !this->analyzeVoxel(vTmp, nl, NULL, bc)){
						newList[index] = vTmp;
						index++;
					}else{
						delete vTmp;
					}

				}while(increment(inangle, ANGULAR_DIMENSION, (unsigned char)1));
			}
		}while(increment(inpos, SPATIAL_DIMENSION, (unsigned char)1));
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

#endif


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> * VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getRandomVoxel(RND *rnd){
	double d = rnd->nextValue();
	return this->voxels[(int)(d*(this->last+1))];
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
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>* VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::get(int i){
	return this->voxels[i];
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::length(){
	return this->last+1;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVoxelsSurface(){
	return (this->last+1)*pow(this->voxelSize, SPATIAL_DIMENSION)*pow((this->angularSize/(2*M_PI)), ANGULAR_DIMENSION);
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toPovray(){
	std::string sRes = "";

	for(int i=0; i<=this->last; i++){
		double *da = this->voxels[i]->getPosition();
		double x1 = da[0], x2 = da[0] + this->voxelSize;
		double y1 = da[1], y2 = da[1] + this->voxelSize;

		sRes += "  polygon {5, < " + std::to_string(x1) + ", " + std::to_string(y1) + ", 1.0>, "
				+ "< " + std::to_string(x1) + ", " + std::to_string(y2) + ", 1.0>, "
				+ "< " + std::to_string(x2) + ", " + std::to_string(y2) + ", 1.0>, "
				+ "< " + std::to_string(x2) + ", " + std::to_string(y1) + ", 1.0>, "
				+ "< " + std::to_string(x1) + ", " + std::to_string(y1) + ", 1.0>"
				+ " texture { pigment { color Black} } }\r\n";
/*
		sRes += "  text { ttf \"timrom.ttf\" \"" + std::to_string(i) + "\" 1, 0 pigment { color White } scale 0.1 translate < ";
		for(unsigned char j=0; j<this->dimension; j++)
			sRes += std::to_string(da[j]+0.5*this->voxelSize) + ", ";
		sRes +=  "1.0003> }\r\n";
*/
	}
	return sRes;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toWolfram(){
	std::stringstream out;

	for(int i=0; i<=this->last; i++){
		double *da = this->voxels[i]->getPosition();
		double x1 = da[0], x2 = da[0] + this->voxelSize;
		double y1 = da[1], y2 = da[1] + this->voxelSize;

		out << "Polygon[{ {" << x1 << ", " << y1 << "}, "
				<< "{" << x1 << ", " << y2 << "}, "
				<< "{" << x2 << ", " << y2 << "}, "
				<< "{" << x2 << ", " << y1 << "} }]";
		if (i!=this->last)
			out << ", ";
		out << std::endl;
		if (ANGULAR_DIMENSION > 0){
			out << "(* angles: [ " << this->voxels[i]->getOrientation()[0] << ", " << (this->voxels[i]->getOrientation()[0] + this->angularSize) << ") *)" << std::endl;
		}
	}
	return out.str();
}
