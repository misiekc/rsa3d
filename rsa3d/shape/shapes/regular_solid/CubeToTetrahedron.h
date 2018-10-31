//
// Created by PKua on 26.09.18.
//

#ifndef RSA3D_CUBETOTETRAHEDRON_H
#define RSA3D_CUBETOTETRAHEDRON_H


#include "RegularSolid.h"
#include "UnoptimizedSATOverlapRS.h"

class CubeToTetrahedron : public RegularSolid<CubeToTetrahedron> {
private:
    friend RegularSolid<CubeToTetrahedron>;

    const static UnoptimizedSATOverlapRS overlapStrategy;

    static void calculateStatic(const std::string &attr);

public:
    explicit CubeToTetrahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;

    static void calculateTruncatedTetrahedron(double xi);

    static void calculateTruncatedCube(double aa);

    static void calculateCuboctahedron(double aa);

    static void calculateNonintersectingTruncations(double aa, double ac);

    static void calculateIntersectingTruncations(double aa, double ac);
};


#endif //RSA3D_CUBETOTETRAHEDRON_H
