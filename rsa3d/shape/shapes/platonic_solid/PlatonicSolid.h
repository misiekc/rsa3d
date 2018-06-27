//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_PLATONICSOLID_H
#define RSA3D_PLATONICSOLID_H


#include "../../../Matrix.h"
#include "../../OverlapStrategyShape.h"
#include "../../../Intersection.h"
#include "../../ConvexShape.h"

// CRTP idiom
template <typename SpecificSolid>
class PlatonicSolid : public ConvexShape<3, 0>, public OverlapStrategyShape<3, 0> {
private:
    Matrix<3, 3> orientation = Matrix<3, 3>::identity();

    static void normalizeAxes();

protected:
    explicit PlatonicSolid(const Matrix<3, 3> &orientation) : orientation(orientation) {};

    void setPosition(const Vector<3> &position) override;

    template <std::size_t SIZE>
    std::array<Vector<3>, SIZE> applyOrientation(const std::array<Vector<3>, SIZE> &vectors) const;
    void calculateVerticesAndAxes();
    void setOrientationMatrix(const Matrix<3, 3> &orientation);
    void translateVertices(const Vector<3> &translation);

public:
    static void initClass(const std::string &attr);

    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const final;
    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                    double orientationRange) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    std::string toWolfram() const override;

    std::vector<std::string> getSupportedStrategies() const override;

    Shape<3, 0> *clone() const override;

    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;

    /* std::array<Vector<3>, ?> getVertices() const;                CRTP pure virtual */
    /* std::array<Vector<3>, ?> getFaceAxes() const;                CRTP pure virtual */
    /* std::array<Vector<3>, ?> getEdgeAxes() const;                CRTP pure virtual */
    /* double projectionHalfsize(const Vector<3> &axis) const;      CRTP pure virtual */
    /* intersection::polyhedron getTriangles() const;               CRTP pure virtual */

    bool isSeparatingAxis(const Vector<3> &axis, const SpecificSolid &other,
                          const Vector<3> &distance) const; /* CRTP virtual */
    const Matrix<3, 3> &getOrientationMatrix() const;
};

#include "PlatonicSolid.tpp"

#endif //RSA3D_PLATONICSOLID_H
