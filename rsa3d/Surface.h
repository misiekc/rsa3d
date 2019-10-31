/*
 * Surface.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef SURFACE_H_
#define SURFACE_H_


#include <vector>

#include "NeighbourGrid.h"
//#include "DeviceNeighbourGrid.cu"

#include "shape/Shape.h"

class Surface : public BoundaryConditions<RSA_SPATIAL_DIMENSION> {

protected:

	NeighbourGrid<const RSAShape> *list;

	double size;
	double dimension;

	RSAVector vectorFreeBC(const RSAVector &v) const;
	RSAVector vectorPeriodicBC(const RSAVector &v) const;

public:
	Surface(int dim, double s, double ndx, double vdx);
	virtual ~Surface();

	void clear();
	void add(const RSAShape *s);
	const RSAShape *check(const RSAShape *s);
	void getNeighbours(std::vector<const RSAShape*> *result, const RSAVector &da);
	const RSAShape *getClosestNeighbour(const RSAVector &da);
	const RSAShape *getClosestNeighbour(const RSAVector &da, const std::vector<const RSAShape*> &neighbours);
	NeighbourGrid<const RSAShape> *getNeighbourGrid();
    double distance2(const RSAVector &a1, const RSAVector &a2) const override;

	RSAVector getTranslation(const RSAVector &p1, const RSAVector &p2) const override = 0;
	virtual RSAVector vector(const RSAVector &v) const = 0;
	virtual RSAVector checkPosition(const RSAVector &da) const;
	virtual double getArea() const;

//	void drawShapes(Graphics g, double scale);
//	void drawShapes(Graphics g, double scale, double[] ta);
};

#endif /* SURFACE_H_ */
