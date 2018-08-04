//
// Created by PKua on 12.07.18.
//

#ifndef RSA3D_TRUNCATEDCUBE_H
#define RSA3D_TRUNCATEDCUBE_H


#include "RegularSolid.h"
#include "../../OrderCalculable.h"

class TruncatedCube : public RegularSolid<TruncatedCube> {
private:
    friend RegularSolid<TruncatedCube>;

    static void calculateStatic(const std::string &attr);

public:
    explicit TruncatedCube(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     double orientationRange) const override;

    double projectionHalfsize(const Vector<3> &axis) const override;
};


#endif //RSA3D_TRUNCATEDCUBE_H
