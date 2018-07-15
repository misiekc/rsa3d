//
// Created by PKua on 15.07.18.
//

#ifndef RSA3D_SNUBDODECAHEDRON_H
#define RSA3D_SNUBDODECAHEDRON_H


#include "RegularSolid.h"

class SnubDodecahedron : public RegularSolid<SnubDodecahedron> {
private:
    friend RegularSolid<SnubDodecahedron>;

    const static TriTriOverlap<SnubDodecahedron> overlapStrategy;

    static void calculateStatic(const std::string &attr);

public:
    explicit SnubDodecahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const;         /* CRTP implement */

    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;
};


#endif //RSA3D_SNUBDODECAHEDRON_H
