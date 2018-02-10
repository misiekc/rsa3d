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

    double angle;

    static double pointDistance2(const Vector<2> &pos, double angle, const Vector<2> &point);
    static double getAngleToOrigin(const Vector<2> &point);

    double pointDistance2(const Vector<2> &p) const;
    bool withinExclusionZone(const Vector<2> &pointPos, double angle);
    Matrix<2, 2> getRotationMatrix() const;
    Matrix<2, 2> getAntiRotationMatrix() const;

public:
    explicit SpheroCylinder2D(double angle);
    SpheroCylinder2D(const SpheroCylinder2D & other) = default;
    ~SpheroCylinder2D() override = default;

    static void initClass(const std::string & attr);
    static Shape<2> * create(RND * rnd);

    double getNeighbourListCellSize() override;
    double getVoxelSize() override;
    int overlap(BoundaryConditions *bc, Shape *s) override;
    double getVolume() override;
    int pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) override;
    int pointInside(BoundaryConditions *bc, double *da) override;
    std::string toString() override;
    std::string toPovray() const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;

    std::string toWolfram() const;
};


#endif //RSA3D_SPHEROCYLINDER2D_H
