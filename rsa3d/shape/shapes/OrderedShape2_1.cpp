//
// Created by pkua on 03.09.2019.
//

#include <sstream>

#include "OrderedShape2_1.h"

/* If defined, zero variation in angle would be actually very small variation - some truly oriented shapes generate
 * very slowly */
#define ORDERED_SHAPE_ZERO_DISPERSION_WORKAROUND

double OrderedShape2_1::preferredAngle{};
double OrderedShape2_1::angleDistributionSigma{};
std::normal_distribution<double> OrderedShape2_1::angleDistribution{};

/* Arguments in format: (angle distribution average) (angle distribution sigma) (encapsulated class parameters) */
void OrderedShape2_1::initClass(const std::string &args, InitClassFunction baseShapeInitClassFunction) {
    std::istringstream argsStream(args);
    argsStream >> preferredAngle >> angleDistributionSigma;
    ValidateMsg(argsStream, "(preferred angle in degrees) (angle sigma in degrees)");
    ValidateMsg(angleDistributionSigma >= 0, "Angle distribution sigma should be non-negative");

#ifdef ORDERED_SHAPE_ZERO_DISPERSION_WORKAROUND
    if (angleDistributionSigma == 0)
        angleDistributionSigma = 0.00000001;
#endif

    preferredAngle = preferredAngle * M_PI / 180.;
    angleDistributionSigma = angleDistributionSigma * M_PI / 180.;
    angleDistribution = std::normal_distribution<double>(preferredAngle, angleDistributionSigma);

    std::string baseShapeArgs;
    std::getline(argsStream, baseShapeArgs);
    baseShapeInitClassFunction(baseShapeArgs);
    auto baseShapeInfo = Shape<2, 1>::getShapeStaticInfo();

    ShapeStaticInfo<2, 0> shapeInfo;
    shapeInfo.setCircumsphereRadius(baseShapeInfo.getCircumsphereRadius());
    shapeInfo.setInsphereRadius(baseShapeInfo.getInsphereRadius());
    shapeInfo.setCreateShapeImpl(createShape);
    Shape<2, 0>::setShapeStaticInfo(shapeInfo);
}

Shape<2, 0> *OrderedShape2_1::createShape(RND *rnd) {
    auto createShape = Shape<2, 1>::getCreateShapeImpl();
    auto underlyingShape = std::unique_ptr<Shape<2, 1>>(createShape(rnd));

    if (angleDistributionSigma > 0) {
        RNDUniformRandomBitGenerator bitGenerator(*rnd);
        double angle = angleDistribution(bitGenerator);
        underlyingShape->rotate({{angle}});
    }
    return new OrderedShape2_1(std::move(underlyingShape));
}

OrderedShape2_1::OrderedShape2_1(std::unique_ptr<Shape<2, 1>> underlyingShape) {
    this->underlyingShape = std::move(underlyingShape);
}

OrderedShape2_1::OrderedShape2_1(const OrderedShape2_1 &other) : Shape(other) {
    this->underlyingShape.reset(other.underlyingShape->clone());
}

OrderedShape2_1 &OrderedShape2_1::operator=(const OrderedShape2_1 &other) {
    Shape::operator=(other);
    this->underlyingShape.reset(other.underlyingShape->clone());
    return *this;
}

bool OrderedShape2_1::overlap(BoundaryConditions<2> *bc, const Shape<2, 0> *s) const {
    auto &orderedOther = dynamic_cast<const OrderedShape2_1&>(*s);
    return this->underlyingShape->overlap(bc, orderedOther.underlyingShape.get());
}

double OrderedShape2_1::getVolume(unsigned short dim) const {
    return this->underlyingShape->getVolume(dim);
}

bool OrderedShape2_1::voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition,
                                  const Orientation<0> &voxelOrientation, double spatialSize,
                                  double angularSize) const {
    return this->underlyingShape->voxelInside(bc, voxelPosition, {{0}}, spatialSize, 2*M_PI);
}

double OrderedShape2_1::minDistance(const Shape<2, 0> *s) const {
    auto &orderedOther = dynamic_cast<const OrderedShape2_1&>(*s);
    return this->underlyingShape->minDistance(orderedOther.underlyingShape.get());
}

std::string OrderedShape2_1::toPovray() const {
    return this->underlyingShape->toPovray();
}

std::string OrderedShape2_1::toWolfram() const {
    return this->underlyingShape->toWolfram();
}

void OrderedShape2_1::store(std::ostream &f) const {
    return this->underlyingShape->store(f);
}

void OrderedShape2_1::restore(std::istream &f) {
    return this->underlyingShape->restore(f);
}

Shape<2, 0> *OrderedShape2_1::clone() const {
    return new OrderedShape2_1(*this);
}

std::string OrderedShape2_1::toString() const {
    return this->underlyingShape->toString();
}

const Vector<2> &OrderedShape2_1::getPosition() const {
    return this->underlyingShape->getPosition();
}

void OrderedShape2_1::translate(const Vector<2> &v) {
    this->underlyingShape->translate(v);
}
