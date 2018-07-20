//
// Created by PKua on 12.07.18.
//

#ifndef RSA3D_TRUNCATEDCUBOCTAHEDRON_H
#define RSA3D_TRUNCATEDCUBOCTAHEDRON_H


#include "RegularSolid.h"

class TruncatedCuboctahedron : public RegularSolid<TruncatedCuboctahedron> {
private:
    friend RegularSolid<TruncatedCuboctahedron>;

    static void calculateStatic(const std::string &attr);

public:
    explicit TruncatedCuboctahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const override;
};


#endif //RSA3D_TRUNCATEDCUBOCTAHEDRON_H
