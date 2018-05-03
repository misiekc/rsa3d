//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARDODECAHEDRON_H
#define RSA3D_REGULARDODECAHEDRON_H


#include "PlatonicSolid.h"

class RegularDodecahedron : public PlatonicSolid<RegularDodecahedron> {
private:
    friend PlatonicSolid<RegularDodecahedron>;

    constexpr static double gold = (1 + std::sqrt(5.)) / 2;
    constexpr static double edge = std::pow(3.75 + 1.75 * std::sqrt(5.), -1./3);
    constexpr static double edgeFactor = edge * gold / 2;

    constexpr static double circumsphereRadius = edge * std::sqrt(3.) / 2 * gold;
    constexpr static double insphereRadius = edge * gold * gold / 2 / std::sqrt(3 - gold);

    static std::array<Vector<3>, 20> orientedVertices;
    static std::array<Vector<3>, 6> orientedFaceAxes;
    static std::array<Vector<3>, 15> orientedEdgeAxes;

    std::array<Vector<3>, 20> vertices;
    std::array<Vector<3>, 6> faceAxes;
    std::array<Vector<3>, 15> edgeAxes;

    static void calculateStatic(const std::string &attr);

public:
    explicit RegularDodecahedron(const Matrix<3, 3> &orientation);

    std::array<Vector<3>, 20> getVertices() const;              /* CRTP implement */
    std::array<Vector<3>, 6> getFaceAxes() const;               /* CRTP implement */
    std::array<Vector<3>, 15> getEdgeAxes() const;              /* CRTP implement */

    double projectionHalfsize(const Vector<3> &axis) const;     /* CRTP implement */
    intersection::polyhedron getTriangles() const;              /* CRTP implement */
};


#endif //RSA3D_REGULARDODECAHEDRON_H
