/*
 * NeighbourGrid.h
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#ifndef NEIGHBOURGRID_H_
#define NEIGHBOURGRID_H_

#include <vector>
#include <set>

#include "Positioned.h"
#include "BoundaryConditions.h"
#include <vector>
#include <unordered_set>

class NeighbourGrid {

private:
	unsigned char dimension;
	double linearSize;
	int n;
	double dx;

	// contains vectors of cells (vectors with Positioned* inside)
	std::vector<std::vector<Positioned* > * > lists;
	std::vector<std::vector<int> * > neighbouringCells;

	void init(int dim, double size, int n);

public:




	NeighbourGrid(int dim, double size, double dx);
	NeighbourGrid(int dim, double size, int n);
	virtual ~NeighbourGrid();

	void add(Positioned* s);
	void remove(Positioned* s);

//	std::unordered_set<Positioned*> * getNeighbours(double* da, unsigned char radius);
	void clear();
	std::vector<Positioned*> * getCell(double* da);
	void getNeighbours(std::unordered_set<Positioned*> *result, double* da);
	Positioned* getClosestNeighbour(double *da, BoundaryConditions *bc);
};

#endif /* NEIGHBOURGRID_H_ */
