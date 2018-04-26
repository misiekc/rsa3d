//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARTETRAHEDRON_H
#define RSA3D_REGULARTETRAHEDRON_H


#include "PlatonicSolid.h"

class RegularTetrahedron : public PlatonicSolid<RegularTetrahedron> {
private:
    friend PlatonicSolid<RegularTetrahedron>;

    constexpr static double circumsphereRadius = std::pow(3, 5./6) / 2;
    constexpr static double insphereRadius = std::pow(3, -1./6) / 2;
    static std::array<Vector<3>, 4> orientedVertices;
    static std::array<Vector<3>, 4> orientedFaceAxes;
    static std::array<Vector<3>, 6> orientedEdgeAxes;

    std::array<Vector<3>, 4> vertices;
    std::array<Vector<3>, 4> faceAxes;
    std::array<Vector<3>, 6> edgeAxes;

    static void calculateStatic(const std::string &attr);

public:
    using interval = std::pair<double, double>;

    explicit RegularTetrahedron(const Matrix<3, 3> &orientation);

    double projectionHalfsize(const Vector<3> &axis) const;
    bool isSeparatingAxis(const Vector<3> &axis, const RegularTetrahedron &other, const Vector<3> &distance) const;
    interval getProjection(const Vector<3> & axis) const;
    std::array<Vector<3>, 4> getVertices() const;
    std::array<Vector<3>, 4> getFaceAxes() const;
    std::array<Vector<3>, 6> getEdgeAxes() const;
};


#endif //RSA3D_REGULARTETRAHEDRON_H
