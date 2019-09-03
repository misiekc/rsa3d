//
// Created by pkua on 03.09.2019.
//

#include "OrderedShape2_1.h"

template<typename BaseShape>
double OrderedShape2_1<BaseShape>::preferredAngle{};

template<typename BaseShape>
double OrderedShape2_1<BaseShape>::angleDistributionSigma{};

template<typename BaseShape>
void OrderedShape2_1<BaseShape>::initClass(const std::string &args) {

}

template<typename BaseShape>
OrderedShape2_1<BaseShape>::OrderedShape2_1(double angle) {

}

template<typename BaseShape>
bool OrderedShape2_1<BaseShape>::overlap(BoundaryConditions<2> *bc, const Shape<2, 0> *s) const {
    return false;
}

template<typename BaseShape>
double OrderedShape2_1<BaseShape>::getVolume(unsigned short dim) const {
    return Shape::getVolume(dim);
}

template<typename BaseShape>
bool OrderedShape2_1<BaseShape>::voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition,
                                             const Orientation<0> &voxelOrientation, double spatialSize,
                                             double angularSize) const {
    return false;
}

template<typename BaseShape>
double OrderedShape2_1<BaseShape>::minDistance(const Shape<2, 0> *s) const {
    return Shape::minDistance(s);
}

template<typename BaseShape>
std::string OrderedShape2_1<BaseShape>::toPovray() const {
    return Shape::toPovray();
}

template<typename BaseShape>
std::string OrderedShape2_1<BaseShape>::toWolfram() const {
    return Shape::toWolfram();
}

template<typename BaseShape>
void OrderedShape2_1<BaseShape>::store(std::ostream &f) const {
    Shape::store(f);
}

template<typename BaseShape>
void OrderedShape2_1<BaseShape>::restore(std::istream &f) {
    Shape::restore(f);
}

template<typename BaseShape>
Shape<2, 0> *OrderedShape2_1<BaseShape>::clone() const {
    return nullptr;
}
