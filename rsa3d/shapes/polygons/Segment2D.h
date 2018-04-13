/*
 * Segment.h
 *
 *  Created on: 29-10-2013
 *      Author: ciesla
 */

#ifndef SEGMENT2D_H_
#define SEGMENT2D_H_

#include "../../Positioned.h"
#include "../../BoundaryConditions.h"

class Segment2D: public Positioned<2>{

private:
	double radius;
	double orientation;
	double x1, x2, y1, y2;

	void calculateXY();

public:
	static double det(double x1, double y1, double x2, double y2, double x3, double y3);
	static int onSegment(double x1, double y1, double x2, double y2, double x3, double y3);
	static int isCrossing(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

	Segment2D();
	Segment2D(double x1, double y1, double x2, double y2);
	int isCrossing(Segment2D *s);
	void applyBC(BoundaryConditions *bc, Segment2D *other) const;
	void translate(double *v) override;
	void rotate(double phi);
	double length();
	double distance(double x, double y);
	double getX1();
	double getY1();
	double getX2();
	double getY2();
};

#endif /* SEGMENT2D_H_ */
