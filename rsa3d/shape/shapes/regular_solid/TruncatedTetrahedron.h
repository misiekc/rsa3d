//
// Created by PKua on 15.07.18.
//

#ifndef RSA3D_TRUNCATEDTETRAHEDRON_H
#define RSA3D_TRUNCATEDTETRAHEDRON_H


#include "RegularSolid.h"

class TruncatedTetrahedron : public RegularSolid<TruncatedTetrahedron> {
private:
    friend RegularSolid<TruncatedTetrahedron>;

    const static TriTriOverlap<TruncatedTetrahedron> overlapStrategy;

    static void calculateStatic(const std::string &attr);

public:
    explicit TruncatedTetrahedron(const Matrix<3, 3> &orientation) : RegularSolid<TruncatedTetrahedron>{orientation} {};

    double projectionHalfsize(const Vector<3> &axis) const;     /* CRTP implement */

    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;
    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const override;

};


#endif //RSA3D_TRUNCATEDTETRAHEDRON_H
