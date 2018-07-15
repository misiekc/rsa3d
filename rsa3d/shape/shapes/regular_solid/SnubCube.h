//
// Created by PKua on 12.07.18.
//

#ifndef RSA3D_SNUBCUBE_H
#define RSA3D_SNUBCUBE_H


#include "RegularSolid.h"
#include "UnoptimizedSATOverlap.h"

class SnubCube : public RegularSolid<SnubCube> {
private:
    friend RegularSolid<SnubCube>;

    static constexpr double tribonacciConstant = (1 + std::cbrt(19 + 3*std::sqrt(33)) + std::cbrt(19 - 3*std::sqrt(33)))/3;
    const static UnoptimizedSATOverlap<SnubCube> overlapStrategy;

    static void calculateStatic(const std::string &attr);

public:
    explicit SnubCube(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    double projectionHalfsize(const Vector<3> &axis) const;         /* CRTP implement */

    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;
};


#endif //RSA3D_SNUBCUBE_H