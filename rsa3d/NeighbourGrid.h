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
#include <vector>
#include <unordered_set>

class NeighbourGrid {

private:
	int* in;

public:

	int dimension;
	double linearSize;
	int n;
	double dx;
	// contains vectors of cells (vectors with Positioned* inside)
	std::vector<std::vector<Positioned* > * > lists;

	std::unordered_set<Positioned *> neighbours;



	NeighbourGrid(int dim, double size, double dx);
	virtual ~NeighbourGrid();

	void add(Positioned* s);
	void remove(Positioned* s);

	std::unordered_set<Positioned*> * getNeighbours(double* da, int radius);
	void clear();
	std::unordered_set<Positioned*> * getNeighbours(double* da);
};

#endif /* NEIGHBOURGRID_H_ */
