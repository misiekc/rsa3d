/*
 * NeighbourGrid.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include "NeighbourGrid.h"
#include <cmath>
#include <algorithm>

#include "Utils.h"

/**
 *
 * @param dim dimension
 * @param size linear size of the structure
 * @param dx linear size of one cell
 */
NeighbourGrid<E>::NeighbourGrid(int dim, double size, double dx) {
	this->init(dim, size, (int)(size/dx));
}

NeighbourGrid<E>::NeighbourGrid(int dim, double size, int n) {
	this->init(dim, size, n);
}

/**
 *
 * @param dim dimension
 * @param size linear size of the structure
 * @param dx linear size of one cell
 */
void NeighbourGrid<E>::init(int dim, double size, int n) {
	this->dimension = dim;
	this->linearSize = size;
	this->n = n;
	this->dx = size/this->n;
	int length = (int) round(pow(this->n, dim));
	this->lists.reserve(length);
	this->neighbouringCells.reserve(length);

	int *in = new int[this->dimension];
	double *da = new double[this->dimension];
	int *coords = new int[this->dimension];

	for(int i=0; i<length; i++){
		this->lists.push_back(new std::vector<Positioned *>);
		this->neighbouringCells.push_back(new std::vector<int>);
		this->neighbouringCells[i]->reserve((1 << this->dimension));

		i2position(da, this->dimension, i, this->dx, this->n);
		for(unsigned char i=0; i<this->dimension; i++){
			in[i] = 0;
		}
		coordinates(coords, da, this->dimension, this->linearSize, this->dx, this->n);
		do{
			int iCell = neighbour2i(coords, in, this->dimension, 1, this->n);
			this->neighbouringCells[i]->push_back(iCell);
		}while(increment(in, this->dimension, 2));
	}
	delete[] coords;
	delete[] da;
	delete[] in;
}

NeighbourGrid::~NeighbourGrid() {
	for(unsigned int i=0; i<this->lists.size(); i++){
		delete this->lists[i];
		delete this->neighbouringCells[i];
	}
}

void NeighbourGrid::add(Positioned* s, double *da){
	int i = position2i(da, this->dimension, this->linearSize, this->dx, this->n);
	this->lists[i]->push_back(s);
}

/*

	protected static <T> boolean compareLists(ArrayList<T> l1, ArrayList<T> l2){
		if(l1.size()!=l2.size())
			return false;
		for(T t: l1){
			if (!l2.contains(t))
				return false;
		}
		return true;
	}
*/

void NeighbourGrid::remove(Positioned* s, double *da){
	int i = position2i(da, this->dimension, this->linearSize, this->dx, this->n);
	std::vector<Positioned *>::iterator it;
	if ( (it = std::find(this->lists[i]->begin(), this->lists[i]->end(), s)) != this->lists[i]->end())
		this->lists[i]->erase(it);
}

std::vector<Positioned*> * NeighbourGrid::getCell(double* da){
	int i = position2i(da, this->dimension, this->linearSize, this->dx, this->n);
	return this->lists[i];
}


void NeighbourGrid::getNeighbours(std::unordered_set<Positioned*> *result, double* da){
	result->clear();
	std::vector<Positioned *> *vTmp;

	int i = position2i(da, this->dimension, this->linearSize, this->dx, this->n);
	for(int iCell : *(this->neighbouringCells[i])){
		vTmp = (this->lists[iCell]);
		result->insert(vTmp->begin(), vTmp->end());
	}
}

/*
std::unordered_set<Positioned*> * NeighbourGrid::getNeighbours(double* da, unsigned char radius){
	this->neighbours.clear();
	std::vector<Positioned *> *vTmp;

	int *in = new int[this->dimension];
	for(unsigned char i=0; i<this->dimension; i++){
		in[i] = 0;
	}

	int *coords = new int[this->dimension];

	coordinates(coords, da, this->dimension, this->linearSize, this->dx, this->n);
	do{
		int i = neighbour2i(coords, in, this->dimension, radius, this->n);
		vTmp = (this->lists[i]);
		this->neighbours.insert(vTmp->begin(), vTmp->end());
	}while(increment(in, this->dimension, 2*radius));
	delete[] coords;
	delete[] in;
	return &this->neighbours;
}
*/

Positioned* NeighbourGrid::getClosestNeighbour(double *da, BoundaryConditions *bc){
	std::vector<Positioned *> *vTmp;

	double d, dmin = std::numeric_limits<double>::max();
	Positioned *pmin = NULL;
	int i = position2i(da, this->dimension, this->linearSize, this->dx, this->n);
	for(int iCell : *(this->neighbouringCells[i])){
		vTmp = (this->lists[iCell]);
		for(Positioned *p : *vTmp){
			d = bc->distance2(da, p->getPosition());
			if (d<dmin){
				pmin = p;
				dmin = d;
			}
		}
	}
	return pmin;
}

void NeighbourGrid::clear(){
	for(unsigned int i=0; i<this->lists.size(); i++){
		this->lists[i]->clear();
	}
}

/*
std::unordered_set<Positioned*> * NeighbourGrid::getNeighbours(double* da){
	return this->getNeighbours(da, 1);
}
*/
