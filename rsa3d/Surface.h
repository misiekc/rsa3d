/*
 * Surface.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef SURFACE_H_
#define SURFACE_H_


#include <vector>
#include <memory>

#include "NeighbourGrid.h"
//#include "DeviceNeighbourGrid.cu"

#include "shape/Shape.h"

class Surface {

protected:
	NeighbourGrid<const RSAShape> *list;
	std::unique_ptr<RSABoundaryConditions> bc;

	double size;
	double dimension;

public:
	Surface(int dim, double s, double ndx, double vdx, std::unique_ptr<RSABoundaryConditions> bc);
	virtual ~Surface();
	Surface(const Surface &) = delete;
	Surface &operator=(const Surface &) = delete;

	void clear();
	void add(const RSAShape *s);
	const RSAShape *check(const RSAShape *s);
	void getNeighbours(std::vector<const RSAShape*> *result, const RSAVector &da);
	const RSAShape *getClosestNeighbour(const RSAVector &da);
	const RSAShape *getClosestNeighbour(const RSAVector &da, const std::vector<const RSAShape*> &neighbours);
	NeighbourGrid<const RSAShape> *getNeighbourGrid();
	RSABoundaryConditions *getBC();

	[[nodiscard]] virtual RSAVector checkPosition(const RSAVector &da) const;
	[[nodiscard]] virtual double getArea() const;
    [[nodiscard]] virtual double getRealArea() const;

//	void drawShapes(Graphics g, double scale);
//	void drawShapes(Graphics g, double scale, double[] ta);
};

#endif /* SURFACE_H_ */
