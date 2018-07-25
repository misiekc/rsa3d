//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_TETRAHEDRON_H
#define RSA3D_TETRAHEDRON_H


#include "RegularSolid.h"
#include "UnoptimizedSATOverlapRS.h"
#include "../../OrderCalculable.h"

class Tetrahedron : public RegularSolid<Tetrahedron>, public OrderCalculable {
private:
    friend RegularSolid<Tetrahedron>;

    const static UnoptimizedSATOverlapRS overlapStrategy;

    static void calculateStatic(const std::string &attr);

public:
    explicit Tetrahedron(const Matrix<3, 3> &orientation) : RegularSolid<Tetrahedron>{orientation} {};

    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     double orientationRange) const override;
    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;

};


#endif //RSA3D_TETRAHEDRON_H
