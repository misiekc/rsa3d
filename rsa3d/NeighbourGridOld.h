/*
 * NeighbourGridOld.h
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#ifndef NEIGHBOURGRIDOLD_H_
#define NEIGHBOURGRIDOLD_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "BoundaryConditions.h"
#include "utils/Utils.h"
#include "utils/Assertions.h"


template <typename E>
class NeighbourGridOld{

private:
	int surfaceDimension;
	double linearSize;
	int n;
	double dx;
	bool withNeighbouringCells;

	// contains vectors of cells (vectors with Positioned* inside)
	std::vector<std::vector<E * > * > lists;
	std::vector<std::vector<int> * > neighbouringCells;


	void fillNeigbouringCellsVector(std::vector<int> *vec, int cellNo) const{
		int in[RSA_SPATIAL_DIMENSION];
		double da[RSA_SPATIAL_DIMENSION];
		int coords[RSA_SPATIAL_DIMENSION];

		vec->reserve((1 << this->surfaceDimension));

		i2position(da, this->surfaceDimension, cellNo, this->dx, this->n);
		for(unsigned char i=0; i<RSA_SPATIAL_DIMENSION; i++){
			in[i] = 0;
		}
		coordinates(coords, da, this->surfaceDimension, this->linearSize, this->dx, this->n);

		do{
			int iCell = neighbour2i(coords, in, this->surfaceDimension, 1, this->n);
			vec->push_back(iCell);
		}while(increment(in, this->surfaceDimension, 2));
		// sort and erase to avoid duplicates - important for small packings
		std::sort( vec->begin(), vec->end() );
		vec->erase( std::unique( vec->begin(), vec->end() ), vec->end() );
	}

	void init(double size, size_t n){
		this->linearSize = size;
		this->n = n;
		this->dx = size/this->n;
		int length = (int) round(pow(this->n, this->surfaceDimension));
		this->lists.reserve(length);
		if (this->withNeighbouringCells)
			this->neighbouringCells.reserve(length);
		else
			this->neighbouringCells.reserve(0);

		for(int i=0; i<length; i++){
			this->lists.push_back(new std::vector<E *>);

			if (this->withNeighbouringCells){
				this->neighbouringCells.push_back(new std::vector<int>);
				this->fillNeigbouringCellsVector(this->neighbouringCells[i], i);
			}
		}
	}

public:

	NeighbourGridOld(int surfaceDim, double size, double dx){ // @suppress("Class members should be properly initialized")
		Expects(surfaceDim > 0);
		Expects(surfaceDim <= RSA_SPATIAL_DIMENSION);
		Expects(size > 0);
		Expects(dx > 0);

		this->surfaceDimension = surfaceDim;
		this->withNeighbouringCells = true;
		size_t n = (size_t)(size/dx);
		if (n == 0)
		    throw std::runtime_error("neighbour grid cell too big");
		this->init(size, n);
	}

	NeighbourGridOld(int surfaceDim, double size, size_t n){ // @suppress("Class members should be properly initialized")
		Expects(surfaceDim > 0);
		Expects(surfaceDim <= RSA_SPATIAL_DIMENSION);
		Expects(size > 0);
		Expects(n > 0);

		this->withNeighbouringCells = true;
		this->surfaceDimension = surfaceDim;
		this->init(size, n);
	}

	NeighbourGridOld(int surfaceDim, double size, size_t n, bool neighbouringCells){ // @suppress("Class members should be properly initialized")
		Expects(surfaceDim > 0);
		Expects(surfaceDim <= RSA_SPATIAL_DIMENSION);
		Expects(size > 0);
		Expects(n > 0);

		this->withNeighbouringCells = neighbouringCells;
		this->surfaceDimension = surfaceDim;
		this->init(size, n);
	}

	virtual ~NeighbourGridOld(){
		for(unsigned int i=0; i<this->lists.size(); i++){
				delete this->lists[i];
				if (withNeighbouringCells)
					delete this->neighbouringCells[i];
			}
	}

	void add(E* s, const RSAVector &da){
	    double array[RSA_SPATIAL_DIMENSION];
	    da.copyToArray(array);
		int i = position2i(array, this->surfaceDimension, this->linearSize, this->dx, this->n);
		this->lists[i]->push_back(s);
	}

	void remove(E* s, const RSAVector &da){
	    double array[RSA_SPATIAL_DIMENSION];
	    da.copyToArray(array);
		int i = position2i(array, this->surfaceDimension, this->linearSize, this->dx, this->n);
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
		int i = position2i(array, this->surfaceDimension, this->linearSize, this->dx, this->n);
		return this->lists[i];
	}

	void getNeighbours(std::vector<E *> *result, const RSAVector &da) const{
		result->clear();
		std::vector<E *> *vTmp;

		double array[RSA_SPATIAL_DIMENSION];
		da.copyToArray(array);
		int i = position2i(array, this->surfaceDimension, this->linearSize, this->dx, this->n);

		std::vector<int> *cellsptr;
		if (this->withNeighbouringCells)
			cellsptr = this->neighbouringCells[i];
		else{
			cellsptr = new std::vector<int>();
			this->fillNeigbouringCellsVector(cellsptr, i);
		}
		for(int iCell : *(cellsptr)){
			vTmp = (this->lists[iCell]);
			result->insert(result->end(), vTmp->begin(), vTmp->end());
		}
		if(!this->withNeighbouringCells)
			delete cellsptr;
	}
};

#endif /* NEIGHBOURGRIDOLD_H_ */
