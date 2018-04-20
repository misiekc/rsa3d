/*
 * Polygon.h
 *
 *  Created on: 16.04.2018
 *      Author: ciesla
 */

#ifndef SHAPES_POLYGONS_POLYGON_H_
#define SHAPES_POLYGONS_POLYGON_H_

#include "../../Shape.h"
#include "Segment2D.h"

class Polygon: public Shape<2, 0> {

private:

	Segment2D **segments;
	unsigned short n;


protected:

	static double findVoxelSpatialSize(Polygon *p);
	static double findNeighbourGridSize(Polygon *p);


public:
	Polygon(unsigned short n);
	Polygon(const Polygon &p);
	virtual ~Polygon();

	void applyBC(BoundaryConditions *bc, Shape<2,0> *second) const override;
	int overlap(BoundaryConditions *bc, Shape<2, 0> *s) const override;
	virtual int pointInside(BoundaryConditions *bc, double* position, const std::array<double, 0> &orientation, double orientationRange) const override;

};

#endif /* SHAPES_POLYGONS_POLYGON_H_ */
