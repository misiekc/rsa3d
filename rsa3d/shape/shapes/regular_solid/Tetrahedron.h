//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_TETRAHEDRON_H
#define RSA3D_TETRAHEDRON_H


#include "RegularSolid.h"
#include "UnoptimizedSATOverlapRS.h"

class CubeToTetrahedron;

class Tetrahedron : public RegularSolid<Tetrahedron> {
private:
    friend RegularSolid<Tetrahedron>;
    friend CubeToTetrahedron;

    const static UnoptimizedSATOverlapRS overlapStrategy;

    static ShapeData calculateStatic(const std::string &attr);

public:
    explicit Tetrahedron(const Matrix<3, 3> &orientation) : RegularSolid<Tetrahedron>{orientation} {};

    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;

};


#endif //RSA3D_TETRAHEDRON_H
