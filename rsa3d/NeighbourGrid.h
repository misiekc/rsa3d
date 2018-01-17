/*
 * NeighbourGrid.h
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#ifndef NEIGHBOURGRID_H_
#define NEIGHBOURGRID_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "BoundaryConditions.h"
#include "Utils.h"


template <typename E>
class NeighbourGrid{

private:
	unsigned char dimension;
	double linearSize;
	int n;
	double dx;

	// contains vectors of cells (vectors with Positioned* inside)
	std::vector<std::vector<E * > * > lists;
	std::vector<std::vector<int> * > neighbouringCells;

	void init(int dim, double size, int n){
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
			this->lists.push_back(new std::vector<E *>);
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

public:

	NeighbourGrid(int dim, double size, double dx){
		this->init(dim, size, (int)(size/dx));
	}


	NeighbourGrid(int dim, double size, int n){
		this->init(dim, size, n);
	}

	virtual ~NeighbourGrid(){
		for(unsigned int i=0; i<this->lists.size(); i++){
				delete this->lists[i];
				delete this->neighbouringCells[i];
			}
	}

	void add(E* s, double *da){
		int i = position2i(da, this->dimension, this->linearSize, this->dx, this->n);
		this->lists[i]->push_back(s);
	}

	void remove(E* s, double *da){
		int i = position2i(da, this->dimension, this->linearSize, this->dx, this->n);
		typename std::vector<E *>::iterator it;
		if ( (it = std::find(this->lists[i]->begin(), this->lists[i]->end(), s)) != this->lists[i]->end())
			this->lists[i]->erase(it);
	}

	void clear(){
		for(unsigned int i=0; i<this->lists.size(); i++){
			this->lists[i]->clear();
		}
	}

	std::vector<E*> * getCell(double* da){
		int i = position2i(da, this->dimension, this->linearSize, this->dx, this->n);
		return this->lists[i];
	}

	void getNeighbours(std::vector<E *> *result, double* da){
		result->clear();
		std::vector<E *> *vTmp;

		int i = position2i(da, this->dimension, this->linearSize, this->dx, this->n);
		for(int iCell : *(this->neighbouringCells[i])){
			vTmp = (this->lists[iCell]);
			result->insert(result->end(), vTmp->begin(), vTmp->end());
		}
	}
};

#endif /* NEIGHBOURGRID_H_ */
