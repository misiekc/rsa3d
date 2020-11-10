/*
 * NeighbourGridNew.h
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
#include "utils/Assertions.h"
#include "utils/Utils.h"


template <typename E>
class NeighbourGrid{


private:
	unsigned short surfaceDimension;
	double linearSize;
	unsigned int n;
	double dx;
	size_t cellSize;

	// contains vectors of cells (vectors with Positioned* inside)
	std::vector<E * > **cells;
	std::vector<int> neighbouringCellsOffsets;


	size_t position2i(const RSAVector &position) const{
		int result = 0;
		int ix;
		for(short i=this->surfaceDimension-1; i>=0; i--){
			Expects(position[i] >= 0);
            // TODO: is this conditions necessary? Maybe this->linearSize + this->dx is also ok?
			Expects(position[i] < this->linearSize);
			ix = (int)(position[i]/this->dx) + 1;

			result = this->n*result + ix;
		}
		return result;
	}

	void coordinates(size_t* result, const RSAVector &position) const{
		for(short i=this->surfaceDimension-1; i>=0; i--){
			Expects(position[i] >= 0);
            // TODO: is this conditions necessary? Maybe this->linearSize + this->dx is also ok?
			Expects(position[i] < this->linearSize);

			result[i] = (int)(position[i]/this->dx) + 1;
		}
	}

	void coordinates(size_t* result, size_t cellNo) const{

		for(unsigned short i=0; i < this->surfaceDimension; i++){
			result[i] = cellNo % this->n;
			cellNo /= this->n;
		}
	}

	size_t coordinates2i(const size_t *coords){
		size_t result = 0;
		for(short i=this->surfaceDimension-1; i>=0; i--){
			result = this->n*result + coords[i];
		}
		return result;
	}

	size_t neighbour2i(size_t* coordinates, int *neighbour, int offset) const{
		size_t result = 0;

		for(short i=this->surfaceDimension-1; i>=0; i--){
			unsigned int ix = coordinates[i] + neighbour[i] - offset;
			Expects(ix >= 0);
			Expects(ix < this->n);
			result = this->n*result + ix;
		}
		return result;
	}

	/*
	 * returns true if cellNo is reflection of a real cell due to periodic boundary conditions
	 */
	bool reflectedCell(size_t cellNo){
		size_t coords[RSA_SPATIAL_DIMENSION];
		this->coordinates(coords, cellNo);
		for(unsigned short i=0; i<this->surfaceDimension; i++){
			if (coords[i]==0 || coords[i]==this->n-1)
				return true;
		}
		return false;
	}
	/*
	 * if cellNo is reflection of a real cell due to periodic boundary conditions
	 * the method returns pointer to the vector in the real cell. Otherwise nullptr is returned
	 */
	std::vector<E *> * reflectedCellVector(size_t cellNo){
		if (!reflectedCell(cellNo))
			return nullptr;

		size_t coords[RSA_SPATIAL_DIMENSION];
		this->coordinates(coords, cellNo);
		for(unsigned short i=0; i<this->surfaceDimension; i++){
			if (coords[i]==0)
				coords[i] = this->n-2;
			if (coords[i]==this->n-1)
				coords[i]=1;
		}
		return this->cells[coordinates2i(coords)];
	}

	void fillNeigbouringCellsOffsets() {
		int in[RSA_SPATIAL_DIMENSION];
		RSAVector position;
		size_t coords[RSA_SPATIAL_DIMENSION];

		this->neighbouringCellsOffsets.reserve((1 << this->surfaceDimension));

		// we are taking the cell somewhere in the middle
		for(unsigned char i=0; i<this->surfaceDimension; i++){
			in[i] = 0;
			coords[i] = this->n/2;
		}
		size_t testCellNo = this->coordinates2i(coords);
		do{
			size_t neigbourNo = this->neighbour2i(coords, in, 1);
			this->neighbouringCellsOffsets.push_back(neigbourNo - testCellNo);
		}while(increment(in, this->surfaceDimension, 2));
		// sort and erase to avoid duplicates - important for small packings
		std::sort( this->neighbouringCellsOffsets.begin(), this->neighbouringCellsOffsets.end() );
		this->neighbouringCellsOffsets.erase( std::unique( this->neighbouringCellsOffsets.begin(), this->neighbouringCellsOffsets.end() ), this->neighbouringCellsOffsets.end() );

	}

	void init(double size, size_t n){
		this->linearSize = size;
		this->n = n;
		this->dx = size/(this->n-2);
		this->cellSize = (int) round(pow(this->n, this->surfaceDimension));
		this->cells = new std::vector<E *> *[this->cellSize];

		// filling real cells
		for(size_t i=0; i<this->cellSize; i++){
			if (this->reflectedCell(i))
				continue;
			this->cells[i] = new std::vector<E *>;
		}

		// filing reflected cells
		for(size_t i=0; i<this->cellSize; i++){
			if (this->reflectedCell(i))
				this->cells[i] = this->reflectedCellVector(i);
		}

		// check
		for(size_t i=0; i<this->cellSize; i++){
			Expects(this->cells[i]!=nullptr);
		}
		this->fillNeigbouringCellsOffsets();
	}

public:

	NeighbourGrid(int surfaceDim, double size, double dx){ // @suppress("Class members should be properly initialized")
		Expects(surfaceDim > 0);
		Expects(surfaceDim <= RSA_SPATIAL_DIMENSION);
		Expects(size > 0);
		Expects(dx > 0);

		this->surfaceDimension = surfaceDim;
		// cells on the edges are used by periodic boundary conditions
		size_t n = (size_t)(size/dx) + 2;
		if (n < 3)
		    throw std::runtime_error("neighbour grid cell too big");
		this->init(size, n);
	}

/*
	NeighbourGrid(int surfaceDim, double size, size_t n){ // @suppress("Class members should be properly initialized")
		Expects(surfaceDim > 0);
		Expects(surfaceDim <= RSA_SPATIAL_DIMENSION);
		Expects(size > 0);
		Expects(n > 0);

		this->surfaceDimension = surfaceDim;
		this->init(size, n);
	}

	NeighbourGridNew(int surfaceDim, double size, size_t n, bool neighbouringCells){ // @suppress("Class members should be properly initialized")
		Expects(surfaceDim > 0);
		Expects(surfaceDim <= RSA_SPATIAL_DIMENSION);
		Expects(size > 0);
		Expects(n > 0);

		this->withNeighbouringCells = neighbouringCells;
		this->surfaceDimension = surfaceDim;
		this->init(size, n);
	}
	*/

	virtual ~NeighbourGrid(){
		for(size_t i=0; i<this->cellSize; i++){
			if (!this->reflectedCell(i))
				delete this->cells[i];
		}
		delete[] this->cells;
	}

	void add(E* s, const RSAVector &position){
		int i = this->position2i(position);
		this->cells[i]->push_back(s);
	}

	void remove(E* s, const RSAVector &position){
		int i = position2i(position);
		typename std::vector<E *>::iterator it;
		if ( (it = std::find(this->cells[i]->begin(), this->cells[i]->end(), s)) != this->cells[i]->end())
			this->cells[i]->erase(it);
	}

	void clear(){
		for(size_t i=0; i<this->cellSize; i++){
			this->cells[i]->clear();
		}
	}

	std::vector<E*> * getCell(const RSAVector &position){
		int i = position2i(position);
		return this->cells[i];
	}

	void getNeighbours(std::vector<E *> *result, const RSAVector &position) const{
		result->clear();
		std::vector<E *> *vTmp;

		int cellNo = this->position2i(position);

		for(int iCellOffset : this->neighbouringCellsOffsets){
			vTmp = (this->cells[cellNo + iCellOffset]);
			result->insert(result->end(), vTmp->begin(), vTmp->end());
		}
	}
};

#endif /* NEIGHBOURGRID_H_ */
