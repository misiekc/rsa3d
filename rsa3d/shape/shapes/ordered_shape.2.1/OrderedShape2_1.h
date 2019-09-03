//
// Created by pkua on 03.09.2019.
//

#ifndef RSA3D_ORDEREDSHAPE2_1_H
#define RSA3D_ORDEREDSHAPE2_1_H

#include "../../Shape.h"

template <typename BaseShape>
class OrderedShape2_1 : public Shape<2, 0> {
private:
    static double preferredAngle;
    static double angleDistributionSigma;

public:
    static void initClass(const std::string &args);

    OrderedShape2_1(double angle);
    ~OrderedShape2_1() override = default;

    bool overlap(BoundaryConditions<2> *bc, const Shape<2, 0> *s) const override;
    double getVolume(unsigned short dim) const override;
    bool voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition, const Orientation<0> &voxelOrientation,
                     double spatialSize, double angularSize) const override;
    double minDistance(const Shape<2, 0> *s) const override;

    std::string toPovray() const override ;
    std::string toWolfram() const override;

    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;

    Shape<2, 0> *clone() const override;
};

#include "OrderedShape2_1.tpp"


#endif //RSA3D_ORDEREDSHAPE2_1_H
