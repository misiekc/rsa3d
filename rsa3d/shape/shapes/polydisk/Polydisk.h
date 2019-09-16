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

#include "../../AnisotropicShape2D.h"
#include "../../../geometry/Vector.h"

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

    static std::vector<std::size_t> getIndicesOfLargestDisks();
    static Vector<2> getApproximateMassCentre();
    static Vector<2> findPositionOfLargestDiskClosestToMassCentre();
    static void centerPolydisk(const Vector<2> &newCentre);
    static Vector<2> getStaticDiskPosition(size_t diskIndex);

    Vector<2> getDiskPosition(std::size_t diskIndex) const;

    //like voxelInside, but optimized for full angle angularSize
    bool fullAngleVoxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition, double spatialSize) const;

protected:
	//polar coordinates of all disks
	//assume vertex 0 is linked to vertex 1, vertex 1 is linked to vertex 2, vertex 2 is linked to vertex 3, etc.
	//assume vertex (VertexR.size()-1) is linked to vertex 0
	static std::vector<double> diskCentreR;
    static std::vector<double> diskCentreTheta;
	static std::vector<double> diskR;
	static double area;

public:
    /**
     * @param args contains information about coordinates of each disk of the Polydisk.
     * Data should be separated by only spaces, and should be:
     * number_of_disks [xy|rt] c01 c02 r0 c11 c12 r1 c21 c22 r2 ... area
     * xy means cartesian coordinates, and rt means polar coordinates
     * Example format of coordinates
     * 2 xy -1 0 1 1 0 1 6.28318530718
     * or equivalently
     * 4 rt 1 0 1 1 3.1415927 1 6.28318530718
     * if provided area == 0 then it is calculated using Monte-Carlo method.
     * Shape is automatically rescaled to have unit area.
     */
	static void initClass(const std::string &args);

    static const std::vector<double> &getDiskCentreR() { return diskCentreR; }
    static const std::vector<double> &getDiskCentreTheta() { return diskCentreTheta; }

    static const std::vector<double> &getDiskR() {
        return diskR;
    }

	Shape<2, 1> *clone() const override;
	double getVolume(unsigned short dim) const override;

	bool overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const override;
	bool voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition, const Orientation<1> &voxelOrientation,
					 double spatialSize, double angularSize) const override;
    std::string toPovray() const override;
    std::string toWolfram() const override;
};

#endif /* SHAPES_POLYDISK_H_ */
