//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_PLATONICSOLID_H
#define RSA3D_PLATONICSOLID_H


#include "../../Shape.h"
#include "../../Matrix.h"
#include "../OverlapStrategyShape.h"
#include "../../Intersection.h"

// CRTP idiom
template <typename SpecificSolid>
class PlatonicSolid : public Shape<3, 0>, public OverlapStrategyShape<3, 0> {
private:
    Matrix<3, 3> orientation = Matrix<3, 3>::identity();

protected:
    explicit PlatonicSolid(const Matrix<3, 3> &orientation) : orientation(orientation) {};

    void setPosition(const double *position) override;

    template <std::size_t SIZE>
    std::array<Vector<3>, SIZE> applyOrientation(const std::array<Vector<3>, SIZE> &vectors) const;
    void calculateVerticesAndAxes();
    void setOrientationMatrix(const Matrix<3, 3> &orientation);

public:
    static void initClass(const std::string &attr);

    int overlap(BoundaryConditions *bc, Shape<3, 0> *s) const final;
    int pointInside(BoundaryConditions *bc, double *position, const std::array<double, 0> &orientation,
                    double orientationRange) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;

    std::vector<std::string> getSupportedStrategies() const override;
    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;

    /* double projectionHalfsize(const Vector<3> &axis) const;      CRTP pure virtual */
    /* std::array<Vector<3>, 4> getVertices() const;                CRTP pure virtual */
    /* std::array<Vector<3>, 4> getFaceAxes() const;                CRTP pure virtual */
    /* std::array<Vector<3>, 6> getEdgeAxes() const;                CRTP pure virtual */

    bool isSeparatingAxis(const Vector<3> &axis, const SpecificSolid &other,
                          const Vector<3> &distance) const; /* CRTP virtual */
    intersection::polyhedron getTriangles() const; /* CRTP pure virtual */

    const Matrix<3, 3> &getOrientationMatrix() const;
};

#include "PlatonicSolid.tpp"

#endif //RSA3D_PLATONICSOLID_H