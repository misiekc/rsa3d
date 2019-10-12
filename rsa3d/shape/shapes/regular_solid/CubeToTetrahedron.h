//
// Created by PKua on 26.09.18.
//

#ifndef RSA3D_CUBETOTETRAHEDRON_H
#define RSA3D_CUBETOTETRAHEDRON_H


#include "RegularSolid.h"
#include "UnoptimizedSATOverlapRS.h"
#include "Octahedron.h"

class CubeToTetrahedron : public RegularSolid<CubeToTetrahedron> {
private:
    friend RegularSolid<CubeToTetrahedron>;

    using SymmetryPlatonicSolid = Octahedron;

    const static UnoptimizedSATOverlapRS overlapStrategy;

    static ShapeData calculateStatic(const std::string &attr);

public:
    explicit CubeToTetrahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;

    static ShapeData calculateTruncatedTetrahedron(double xi);
    static ShapeData calculateTruncatedCube(double aa);
    static ShapeData calculateCuboctahedron(double aa);
    static ShapeData calculateNonintersectingTruncations(double aa, double ac);
    static ShapeData calculateIntersectingTruncations(double aa, double ac);
};


#endif //RSA3D_CUBETOTETRAHEDRON_H
