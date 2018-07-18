//
// Created by PKua on 17.07.18.
//

#ifndef RSA3D_ELLIPSOID_H
#define RSA3D_ELLIPSOID_H


#include "../ConvexShape.h"

class Ellipsoid : public ConvexShape<3, 0> {
private:
    static void normalizeVolume();
    static double a, b, c;

    Matrix<3, 3> orientation;

public:
    static void initClass(const std::string &attr);

    explicit Ellipsoid(const Matrix<3, 3> &orientation) : orientation(orientation) {}

    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const override;
    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     double orientationRange) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    Shape<3, 0> *clone() const override;

    std::string toPovray() const;
};


#endif //RSA3D_ELLIPSOID_H
