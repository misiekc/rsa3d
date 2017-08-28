/*
 * Ellipse.h
 *
 *  Created on: 10.08.2017
 *      Author: ciesla
 */

#ifndef SHAPES_ELLIPSE_H_
#define SHAPES_ELLIPSE_H_

#include "../Shape.h"
#include "../RND.h"
#include "../BoundaryConditions.h"

class Ellipse: public Shape {
private:
	static double longSemiAxis;
	static double shortSemiAxis;
	static double neighbourListCellSize;
	static double voxelSize;

	static void rotate(double* point, double alpha);

	double a, b;
	double orientation;
	double u[2], uT[2];

	void calculateU();
	Ellipse & operator = (const Ellipse & el);

	double calculateF(double* r, double g);


public:
	Ellipse();
	virtual ~Ellipse();

	static void initClass(const std::string &args);
	static Shape * create(RND *rnd);

	int overlap(BoundaryConditions *bc, Shape *s);
	double getVolume();
	int pointInside(BoundaryConditions *bc, double* da);
	int pointInside(BoundaryConditions *bc, double* da, double angleFrom, double angleTo);
	double getNeighbourListCellSize();
	double getVoxelSize();

	std::string toPovray() const;
	void store(std::ostream &f) const;
	void restore(std::istream &f);

};

#endif /* SHAPES_ELLIPSE_H_ */