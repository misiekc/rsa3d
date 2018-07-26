//
// Created by PKua on 14.07.18.
//

#ifndef RSA3D_RHOMBICOSIDODECAHEDRON_H
#define RSA3D_RHOMBICOSIDODECAHEDRON_H


#include "RegularSolid.h"
#include "../../OrderCalculable.h"

class Rhombicosidodecahedron : public RegularSolid<Rhombicosidodecahedron> {
private:
    friend RegularSolid<Rhombicosidodecahedron>;

    static void calculateStatic(const std::string &attr);

public:
    explicit Rhombicosidodecahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const override;
};


#endif //RSA3D_RHOMBICOSIDODECAHEDRON_H
