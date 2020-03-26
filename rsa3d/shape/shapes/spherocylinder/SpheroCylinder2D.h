//
// Created by PKua on 07.02.18.
//

#ifndef RSA3D_SPHEROCYLINDER2D_H
#define RSA3D_SPHEROCYLINDER2D_H


#include "../../AnisotropicShape2D.h"
#include "../../../geometry/Vector.h"

class SpheroCylinder2D : public AnisotropicShape2D {
private:
    static double radius;
    static double halfDistance;
    static Vector<2> centerVector;

    /**
     * @brief A ShapeStaticInfo for SpheroCylinder2D. This one is used, instead of Shape class static info, because
     * spherocylinder in also used in RegularDiskoPolygon, which overrides Shape's static info.
     */
    static ShapeStaticInfo<2, 1> spherocylinderShapeInfo;

    static double pointDistance2(const Vector<2> &pos, double angle, const Vector<2> &point);

    bool angleInRange(double angle, double rangeStart, double rangeEnd) const;
    double pointDistance2(const Vector<2> &p) const;
    bool withinExclusionZone(const Vector<2> &pointPos, double angle) const;

protected:
    void setAngle(double angle) override;

public:
    static void initClass(const std::string & attr);
    static void calculateStatic(const std::string &attr);

    double getVolume(unsigned short dim) const override;

    Shape<2, 1> *clone() const override;

    bool overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const override;
    bool pointInside(BoundaryConditions<2> *bc, const Vector<2> &da, double angleFrom, double angleTo) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    std::string toString() const override;
    std::string toPovray() const override;
    std::string toWolfram() const override;
};


#endif //RSA3D_SPHEROCYLINDER2D_H
