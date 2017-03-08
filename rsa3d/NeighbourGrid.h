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

class NeighbourGrid {
public:

	int dimension;
	double linearSize;
	int n;
	double dx;
	std::vector<std::set<Shape*>> lists;



	NeighbourGrid(int dim, double size, double dx);
	virtual ~NeighbourGrid();

	void add(Shape* s);
	void remove(Shape* s);

	std::vector<Shape*> getNeighbours(double* da, int radius);
	void NeighbourGrid::clear();
	std::vector<Shape*> NeighbourGrid::getNeighbours(double* da);
};

#endif /* NEIGHBOURGRID_H_ */
