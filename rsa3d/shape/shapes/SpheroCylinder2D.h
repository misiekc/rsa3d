//
// Created by PKua on 07.02.18.
//

#ifndef RSA3D_SPHEROCYLINDER2D_H
#define RSA3D_SPHEROCYLINDER2D_H


#include "../AnisotropicShape2D.h"
#include "../../Vector.h"

class SpheroCylinder2D : public AnisotropicShape2D {
private:
    static double radius;
    static double halfDistance;
    static Vector<2> centerVector;

    static double pointDistance2(const Vector<2> &pos, double angle, const Vector<2> &point);

    bool angleInRange(double angle, double rangeStart, double rangeEnd) const;
    double pointDistance2(const Vector<2> &p) const;
    bool withinExclusionZone(const Vector<2> &pointPos, double angle) const;

public:
    static void initClass(const std::string & attr);

    double getVolume() const override;

    Shape<2, 1> *clone() const override;

    bool overlap(BoundaryConditions *bc, Shape<2, 1> *s) const override;
    bool pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    std::string toString() const override;
    std::string toPovray() const override;
    std::string toWolfram() const override;
};


#endif //RSA3D_SPHEROCYLINDER2D_H
