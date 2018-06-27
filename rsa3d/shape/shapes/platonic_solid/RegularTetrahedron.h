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

    std::array<Vector<3>, 4> getVertices() const;  /* CRTP implement */
    std::array<Vector<3>, 4> getFaceAxes() const;  /* CRTP implement */
    std::array<Vector<3>, 6> getEdgeAxes() const;  /* CRTP implement */

    double projectionHalfsize(const Vector<3> &axis) const;     /* CRTP implement */
    intersection::polyhedron getTriangles() const;              /* CRTP implement */
    bool isSeparatingAxis(const Vector<3> &axis, const RegularTetrahedron &other,
                          const Vector<3> &distance) const;     /* CRTP override */

    interval getProjection(const Vector<3> & axis) const;

    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     double orientationRange) const override;
};


#endif //RSA3D_REGULARTETRAHEDRON_H
