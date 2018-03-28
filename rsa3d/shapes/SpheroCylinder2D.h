//
// Created by PKua on 07.02.18.
//

#ifndef RSA3D_SPHEROCYLINDER2D_H
#define RSA3D_SPHEROCYLINDER2D_H


#include "../AnisotropicShape2D.h"
#include "../Vector.h"

class SpheroCylinder2D : public AnisotropicShape2D {
private:
    static double voxelSize;
    static double neighbourListCellSize;
    static double radius;
    static double halfDistance;
    static Vector<2> centerVector;

    static double pointDistance2(const Vector<2> &pos, double angle, const Vector<2> &point);

    bool angleInRange(double angle, double rangeStart, double rangeEnd) const;
    double pointDistance2(const Vector<2> &p) const;
    bool withinExclusionZone(const Vector<2> &pointPos, double angle);

public:
    static void initClass(const std::string & attr);
    static Shape<2, 1> * create(RND * rnd);

    double getVoxelAngularSize() override;
    double getNeighbourListCellSize() override;
    double getVoxelSize() override;
    double getVolume() override;
    int overlap(BoundaryConditions *bc, Shape<2, 1> *s) override;
    int pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    std::string toString() override;
    std::string toPovray() const override;
    std::string toWolfram() const override;
};


#endif //RSA3D_SPHEROCYLINDER2D_H
