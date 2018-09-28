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
	double linearSize;
	int n;
	double dx;

	// contains vectors of cells (vectors with Positioned* inside)
	std::vector<std::vector<E * > * > lists;
	std::vector<std::vector<int> * > neighbouringCells;

	void init(double size, size_t n){
		this->linearSize = size;
		this->n = n;
		this->dx = size/this->n;
		int length = (int) round(pow(this->n, RSA_SPATIAL_DIMENSION));
		this->lists.reserve(length);
		this->neighbouringCells.reserve(length);

		int in[RSA_SPATIAL_DIMENSION];
		double da[RSA_SPATIAL_DIMENSION];
		int coords[RSA_SPATIAL_DIMENSION];

		for(int i=0; i<length; i++){
			this->lists.push_back(new std::vector<E *>);
			this->neighbouringCells.push_back(new std::vector<int>);
			this->neighbouringCells[i]->reserve((1 << RSA_SPATIAL_DIMENSION));

			i2position(da, RSA_SPATIAL_DIMENSION, i, this->dx, this->n);
			for(unsigned char i=0; i<RSA_SPATIAL_DIMENSION; i++){
				in[i] = 0;
			}
			coordinates(coords, da, RSA_SPATIAL_DIMENSION, this->linearSize, this->dx, this->n);
			do{
				int iCell = neighbour2i(coords, in, RSA_SPATIAL_DIMENSION, 1, this->n);
				this->neighbouringCells[i]->push_back(iCell);
			}while(increment(in, RSA_SPATIAL_DIMENSION, 2));
			// sort and erase to avoid duplicates - important for small packings
			std::sort( neighbouringCells[i]->begin(), neighbouringCells[i]->end() );
			neighbouringCells[i]->erase( std::unique( neighbouringCells[i]->begin(), neighbouringCells[i]->end() ), neighbouringCells[i]->end() );
		}
	}

public:

	NeighbourGrid(double size, double dx){ // @suppress("Class members should be properly initialized")
		if (size <= 0 || dx <= 0)
		    throw std::runtime_error("size <= 0 || dx <= 0");

		size_t n = (size_t)(size/dx);
		if (n == 0)
		    throw std::runtime_error("neighbour grid cell too big");
		this->init(size, n);
	}


	NeighbourGrid(double size, size_t n){ // @suppress("Class members should be properly initialized")
		this->init(size, n);
	}

	virtual ~NeighbourGrid(){
		for(unsigned int i=0; i<this->lists.size(); i++){
				delete this->lists[i];
				delete this->neighbouringCells[i];
			}
	}

	void add(E* s, const RSAVector &da){
	    double array[RSA_SPATIAL_DIMENSION];
	    da.copyToArray(array);
		int i = position2i(array, RSA_SPATIAL_DIMENSION, this->linearSize, this->dx, this->n);
		this->lists[i]->push_back(s);
	}

	void remove(E* s, const RSAVector &da){
	    double array[RSA_SPATIAL_DIMENSION];
	    da.copyToArray(array);
		int i = position2i(array, RSA_SPATIAL_DIMENSION, this->linearSize, this->dx, this->n);
		typename std::vector<E *>::iterator it;
		if ( (it = std::find(this->lists[i]->begin(), this->lists[i]->end(), s)) != this->lists[i]->end())
			this->lists[i]->erase(it);
	}

	void clear(){
		for(unsigned int i=0; i<this->lists.size(); i++){
			this->lists[i]->clear();
		}
	}

	std::vector<E*> * getCell(const RSAVector &da){
	    double array[RSA_SPATIAL_DIMENSION];
	    da.copyToArray(array);
		int i = position2i(array, RSA_SPATIAL_DIMENSION, this->linearSize, this->dx, this->n);
		return this->lists[i];
	}

	void getNeighbours(std::vector<E *> *result, const RSAVector &da) const{
		result->clear();
		std::vector<E *> *vTmp;

		double array[RSA_SPATIAL_DIMENSION];
		da.copyToArray(array);
		int i = position2i(array, RSA_SPATIAL_DIMENSION, this->linearSize, this->dx, this->n);
		for(int iCell : *(this->neighbouringCells[i])){
			vTmp = (this->lists[iCell]);
			result->insert(result->end(), vTmp->begin(), vTmp->end());
		}
	}
};

#endif /* NEIGHBOURGRID_H_ */
