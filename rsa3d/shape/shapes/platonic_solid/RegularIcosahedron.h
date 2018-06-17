//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARICOSAHEDRON_H
#define RSA3D_REGULARICOSAHEDRON_H


#include "PlatonicSolid.h"

class RegularIcosahedron : public PlatonicSolid<RegularIcosahedron> {
private:
    friend PlatonicSolid<RegularIcosahedron>;

    constexpr static double gold = (1 + std::sqrt(5.)) / 2;
    constexpr static double edge = std::pow(1.25 + 5. / 12 * std::sqrt(5.), -1./3);
    constexpr static double edgeFactor = edge / 2;

    constexpr static double circumsphereRadius = edge * std::sqrt(gold * std::sqrt(5)) / 2;
    constexpr static double insphereRadius = edge * gold * gold / 2 / std::sqrt(3);
    
    static std::array<Vector<3>, 20> orientedVertices;
    static std::array<Vector<3>, 10> orientedFaceAxes;
    static std::array<Vector<3>, 15> orientedEdgeAxes;

    std::array<Vector<3>, 20> vertices;
    std::array<Vector<3>, 10> faceAxes;
    std::array<Vector<3>, 15> edgeAxes;

    static void calculateStatic(const std::string &attr);

public:
    explicit RegularIcosahedron(const Matrix<3, 3> &orientation);

    std::array<Vector<3>, 20> getVertices() const;              /* CRTP implement */
    std::array<Vector<3>, 10> getFaceAxes() const;              /* CRTP implement */
    std::array<Vector<3>, 15> getEdgeAxes() const;              /* CRTP implement */
    double projectionHalfsize(const Vector<3> &axis) const;     /* CRTP implement */
    intersection::polyhedron getTriangles() const;              /* CRTP implement */
};


#endif //RSA3D_REGULARICOSAHEDRON_H
