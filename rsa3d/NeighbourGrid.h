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

class NeighbourGrid {
public:

	int dimension;
	double linearSize;
	int n;
	double dx;
	std::vector<std::set<Positioned*>> lists;



	NeighbourGrid(int dim, double size, double dx);
	virtual ~NeighbourGrid();

	void add(Positioned* s);
	void remove(Positioned* s);

	std::vector<Positioned*> getNeighbours(double* da, int radius);
	void NeighbourGrid::clear();
	std::vector<Positioned*> NeighbourGrid::getNeighbours(double* da);
};

#endif /* NEIGHBOURGRID_H_ */
