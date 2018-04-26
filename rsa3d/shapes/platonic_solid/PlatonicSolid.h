//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_PLATONICSOLID_H
#define RSA3D_PLATONICSOLID_H


#include "../../Shape.h"
#include "../../Matrix.h"
#include "../OverlapStrategyShape.h"

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

    void calculateVertices();
    void calculateAxes();
    void setOrientationMatrix(const Matrix<3, 3> &orientation);
    bool isSeparatingAxis(const Vector<3> &axis, const SpecificSolid &other, const Vector<3> &distance) const;

public:
    static void initClass(const std::string &attr);

    int overlap(BoundaryConditions *bc, Shape<3, 0> *s) const override;
    int pointInside(BoundaryConditions *bc, double *position, const std::array<double, 0> &orientation,
                    double orientationRange) const override;
    void store(std::ostream &f) const override;

    void restore(std::istream &f) override;

    const Matrix<3, 3> &getOrientationMatrix() const;

    std::vector<std::string> getSupportedStrategiesNames() const override;

    OverlapStrategy<3, 0> *getStrategy(const std::string &name) const override;
};

#include "PlatonicSolid.tpp"

#endif //RSA3D_PLATONICSOLID_H
