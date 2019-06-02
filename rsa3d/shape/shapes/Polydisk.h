/*
 * Polydisk.h
 *
 *  Created on: 23.09.2018
 *      Author: ciesla
 */

#ifndef SHAPES_POLYDISK_H_
#define SHAPES_POLYDISK_H_

#include <vector>
#include <utility>
#include <cstddef>

#include "../AnisotropicShape2D.h"
#include "../../geometry/Vector.h"

class Polydisk : public Shape<2, 1> {

private:

	// finds shape area using Monte-Carlo sampling
	static double mcArea(size_t mcTrials);

	// normalize shape to have unit area
	static void normalizeArea(double area);

	// finds minimum and maximum value of cosine function in [theta, theta+dt]
	static Vector<2> minmaxCos(double theta, double dt);

	// finds minimum and maximum value of sine function in [theta, theta+dt]
	static Vector<2> minmaxSin(double theta, double dt);

	//test if line segment from point 1 to 2 intersects with line segment from point 3 to 4
	static bool diskDiskIntersect(size_t disk0, const Vector<2> disk0Position, const Orientation<1> disk0Orientation,
								  size_t disk1, const Vector<2> disk1Position, const Orientation<1> disk1Orientation);

	//same as above, except that endpoints 3 and 4 comes from a line in a voxel, and thus carry an uncertainty
	static bool diskVoxelIntersect(size_t disk0, const Vector<2> disk0Position, const Orientation<1> disk0Orientation,
			 	 	 	 	 	   size_t disk1, const Vector<2> disk1Position, const Orientation<1> disk1Orientation,
								   double spatialSize, double angularSize);

protected:
	//polar coordinates of all disks
	//assume vertex 0 is linked to vertex 1, vertex 1 is linked to vertex 2, vertex 2 is linked to vertex 3, etc.
	//assume vertex (VertexR.size()-1) is linked to vertex 0
	static std::vector<double> diskCentreR;
	static std::vector<double> diskCentreTheta;
	static std::vector<double> diskR;
	static double area;

public:
	static void initClass(const std::string &args);

	~Polydisk() override = default;

	Shape<2, 1> *clone() const override;
	double getVolume(unsigned short dim) const override;

	bool overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const override;
	bool voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition, const Orientation<1> &voxelOrientation,
					 double spatialSize, double angularSize) const override;
	std::string toPovray() const override;
};

#endif /* SHAPES_POLYDISK_H_ */
