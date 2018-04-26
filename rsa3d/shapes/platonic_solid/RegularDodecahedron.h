//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARDODECAHEDRON_H
#define RSA3D_REGULARDODECAHEDRON_H


#include "PlatonicSolid.h"

class RegularDodecahedron : public PlatonicSolid<RegularDodecahedron> {
private:
    friend PlatonicSolid<RegularDodecahedron>;

    constexpr static double circumsphereRadius = 0;
    constexpr static double insphereRadius = 0;
    static std::array<Vector<3>, 12> orientedVertices;
    static std::array<Vector<3>, 6> orientedFaceAxes;
    static std::array<Vector<3>, 15> orientedEdgeAxes;

    std::array<Vector<3>, 12> vertices;
    std::array<Vector<3>, 6> faceAxes;
    std::array<Vector<3>, 15> edgeAxes;

    static void calculateStatic(const std::string &attr);

public:
    explicit RegularDodecahedron(const Matrix<3, 3> &orientation);

    double projectionHalfsize(const Vector<3> &axis) const;
    std::array<Vector<3>, 12> getVertices() const;
    std::array<Vector<3>, 6> getFaceAxes() const;
    std::array<Vector<3>, 15> getEdgeAxes() const;
};


#endif //RSA3D_REGULARDODECAHEDRON_H
