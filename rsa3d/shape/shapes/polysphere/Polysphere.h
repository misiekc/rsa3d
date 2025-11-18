//
// Created by ciesla on 11/13/25.
//

#ifndef RSA3D_POLYSPHERE_H
#define RSA3D_POLYSPHERE_H

#include"../../Shape.h"

class Polysphere : public Shape<3, 3>{
private:

    // finds shape area using Monte-Carlo sampling
    static double mcArea(size_t mcTrials);

    // normalize shape to have unit area
    static void normalizeArea(double area);

    // finds minimum and maximum value of cosine function in [theta, theta+dt]
    static std::array<double, 2> minmaxCos(double theta, double dt);

    // finds minimum and maximum value of sine function in [theta, theta+dt]
    static std::array<double, 2> minmaxSin(double theta, double dt);

    //test if line segment from point 1 to 2 intersects with line segment from point 3 to 4
    static bool sphereSphereIntersect(size_t sphere0, const Vector<3> sphere0Position, const Orientation<3> sphere0Orientation,
                                  size_t sphere1, const Vector<3> sphere1Position, const Orientation<3> sphere1Orientation);

    //same as above, except that endpoints 3 and 4 comes from a line in a voxel, and thus carry an uncertainty
    static bool sphereVoxelIntersect(size_t sphere0, const Vector<3> sphere0Position, const Orientation<3> sphere0Orientation,
                                   size_t sphere1, const Vector<3> sphere1Position, const Orientation<3> sphere1Orientation,
                                   double spatialSize, double angularSize, std::array<std::array<double, 2>, 6> minMaxTrigonometricArray);

    static std::array<double, 3> getStaticSpherePosition(size_t diskIndex, Vector<3> position, Orientation<3> orientation);
    std::array<double, 3> getSpherePosition(size_t index) const;


        //like voxelInside, but optimized for full angle angularSize
    bool fullAngleVoxelInside(BoundaryConditions<3> *bc, const Vector<3> &voxelPosition, double spatialSize) const;

protected:
    //polar coordinates of all disks
    //assume vertex 0 is linked to vertex 1, vertex 1 is linked to vertex 2, vertex 2 is linked to vertex 3, etc.
    //assume vertex (VertexR.size()-1) is linked to vertex 0
    static std::vector<std::array<double, 3>> sphereCentre;
    static std::vector<double> sphereR;

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

    static const std::vector<std::array<double, 3>> &getSphereCentre() { return sphereCentre; }

    static const std::vector<double> &getDiskR() {
        return sphereR;
    }

    Shape<3, 3> *clone() const override;
    double getVolume(unsigned short dim) const override;

    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 3> *s) const override;
    bool voxelInside(BoundaryConditions<3> *bc, const Vector<3> &voxelPosition, const Orientation<3> &voxelOrientation,
                     double spatialSize, double angularSize) const override;
    std::string toPovray() const override;
    std::string toWolfram() const override;
};

#endif //RSA3D_POLYSPHERE_H
