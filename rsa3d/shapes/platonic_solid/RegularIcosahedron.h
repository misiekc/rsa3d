//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARICOSAHEDRON_H
#define RSA3D_REGULARICOSAHEDRON_H


#include "PlatonicSolid.h"

class RegularIcosahedron : public PlatonicSolid<RegularIcosahedron> {
private:
    friend PlatonicSolid<RegularIcosahedron>;

    constexpr static double circumsphereRadius = 0;
    constexpr static double insphereRadius = 0;
    static std::array<Vector<3>, 20> orientedVertices;
    static std::array<Vector<3>, 10> orientedFaceAxes;
    static std::array<Vector<3>, 15> orientedEdgeAxes;

    std::array<Vector<3>, 20> vertices;
    std::array<Vector<3>, 10> faceAxes;
    std::array<Vector<3>, 15> edgeAxes;

    static void calculateStatic(const std::string &attr);

public:
    explicit RegularIcosahedron(const Matrix<3, 3> &orientation);

    double projectionHalfsize(const Vector<3> &axis) const;
    std::array<Vector<3>, 20> getVertices() const;
    std::array<Vector<3>, 10> getFaceAxes() const;
    std::array<Vector<3>, 15> getEdgeAxes() const;
};


#endif //RSA3D_REGULARICOSAHEDRON_H
