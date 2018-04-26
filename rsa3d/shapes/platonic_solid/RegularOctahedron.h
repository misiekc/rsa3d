//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULAROCTAHEDRON_H
#define RSA3D_REGULAROCTAHEDRON_H


#include "PlatonicSolid.h"

class RegularOctahedron : public PlatonicSolid<RegularOctahedron> {
private:
    friend PlatonicSolid<RegularOctahedron>;

    constexpr static double circumsphereRadius = 0;
    constexpr static double insphereRadius = 0;
    static std::array<Vector<3>, 8> orientedVertices;
    static std::array<Vector<3>, 4> orientedFaceAxes;
    static std::array<Vector<3>, 6> orientedEdgeAxes;

    std::array<Vector<3>, 8> vertices;
    std::array<Vector<3>, 4> faceAxes;
    std::array<Vector<3>, 6> edgeAxes;

    static void calculateStatic(const std::string &attr);

public:
    explicit RegularOctahedron(const Matrix<3, 3> &orientation);

    double projectionHalfsize(const Vector<3> &axis) const;
    std::array<Vector<3>, 8> getVertices() const;
    std::array<Vector<3>, 4> getFaceAxes() const;
    std::array<Vector<3>, 6> getEdgeAxes() const;
};


#endif //RSA3D_REGULAROCTAHEDRON_H
