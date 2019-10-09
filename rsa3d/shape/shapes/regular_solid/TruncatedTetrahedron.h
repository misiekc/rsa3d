//
// Created by PKua on 15.07.18.
//

#ifndef RSA3D_TRUNCATEDTETRAHEDRON_H
#define RSA3D_TRUNCATEDTETRAHEDRON_H


#include "RegularSolid.h"
#include "UnoptimizedSATOverlapRS.h"
#include "../../OrderCalculable.h"

class TruncatedTetrahedron : public RegularSolid<TruncatedTetrahedron> {
private:
    friend RegularSolid<TruncatedTetrahedron>;

    const static UnoptimizedSATOverlapRS overlapStrategy;

    static ShapeData calculateStatic(const std::string &attr);

public:
    explicit TruncatedTetrahedron(const Matrix<3, 3> &orientation) : RegularSolid<TruncatedTetrahedron>{orientation} {};

    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;
};


#endif //RSA3D_TRUNCATEDTETRAHEDRON_H
