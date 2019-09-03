//
// Created by pkua on 03.09.2019.
//

#include <sstream>


template<typename BaseShape>
double OrderedShape2_1<BaseShape>::preferredAngle{};

template<typename BaseShape>
double OrderedShape2_1<BaseShape>::angleDistributionSigma{};

template<typename BaseShape>
std::normal_distribution<double> OrderedShape2_1<BaseShape>::angleDistribution{};


template<typename BaseShape>
void OrderedShape2_1<BaseShape>::initClass(const std::string &args) {
    std::istringstream argsStream(args);
    argsStream >> preferredAngle >> angleDistributionSigma;
    ValidateMsg(argsStream, "Malformed OrderedShape2_1 parameters");
    ValidateMsg(angleDistributionSigma > 0, "Angle distribution sigma should be positive");
    angleDistribution = std::normal_distribution<double>(preferredAngle, angleDistributionSigma);

    std::string baseShapeArgs;
    std::getline(argsStream, baseShapeArgs);
    BaseShape::initClass(baseShapeArgs);
    auto baseShapeInfo = Shape<2, 1>::getShapeStaticInfo();

    ShapeStaticInfo<2, 0> shapeInfo;
    shapeInfo.setCircumsphereRadius(baseShapeInfo.getCircumsphereRadius());
    shapeInfo.setInsphereRadius(baseShapeInfo.getInsphereRadius());
    shapeInfo.setCreateShapeImpl(createShape);
    Shape<2, 0>::setShapeStaticInfo(shapeInfo);
}

template<typename BaseShape>
Shape<2, 0> *OrderedShape2_1<BaseShape>::createShape(RND *rnd) {
    auto createShape = BaseShape::getCreateShapeImpl();
    auto underlyingShape = std::unique_ptr<Shape<2, 1>>(createShape(rnd));

    RNDUniformRandomBitGenerator bitGenerator(*rnd);
    double angle = angleDistribution(bitGenerator);

    underlyingShape->rotate({{angle}});
    return new OrderedShape2_1(std::move(underlyingShape));
}

template<typename BaseShape>
OrderedShape2_1<BaseShape>::OrderedShape2_1(std::unique_ptr<Shape<2, 1>> underlyingShape) {
    this->underlyingShape = std::move(underlyingShape);
}

template<typename BaseShape>
bool OrderedShape2_1<BaseShape>::overlap(BoundaryConditions<2> *bc, const Shape<2, 0> *s) const {
    auto &orderedOther = dynamic_cast<const OrderedShape2_1&>(*s);
    return this->underlyingShape->overlap(bc, orderedOther.underlyingShape.get());
}

template<typename BaseShape>
double OrderedShape2_1<BaseShape>::getVolume(unsigned short dim) const {
    return this->underlyingShape->getVolume(dim);
}

template<typename BaseShape>
bool OrderedShape2_1<BaseShape>::voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition,
                                             const Orientation<0> &voxelOrientation, double spatialSize,
                                             double angularSize) const {
    return this->underlyingShape->voxelInside(bc, voxelPosition, {{0}}, spatialSize, 2*M_PI);
}

template<typename BaseShape>
double OrderedShape2_1<BaseShape>::minDistance(const Shape<2, 0> *s) const {
    auto &orderedOther = dynamic_cast<const OrderedShape2_1&>(*s);
    return this->underlyingShape->minDistance(orderedOther.underlyingShape.get());
}

template<typename BaseShape>
std::string OrderedShape2_1<BaseShape>::toPovray() const {
    return this->underlyingShape->toPovray();
}

template<typename BaseShape>
std::string OrderedShape2_1<BaseShape>::toWolfram() const {
    return this->underlyingShape->toWolfram();
}

template<typename BaseShape>
void OrderedShape2_1<BaseShape>::store(std::ostream &f) const {
    return this->underlyingShape->store(f);
}

template<typename BaseShape>
void OrderedShape2_1<BaseShape>::restore(std::istream &f) {
    return this->underlyingShape->restore(f);
}

template<typename BaseShape>
Shape<2, 0> *OrderedShape2_1<BaseShape>::clone() const {
    auto underlyingShapeCopy = this->underlyingShape->clone();
    return new OrderedShape2_1(std::unique_ptr<Shape<2, 1>>(underlyingShapeCopy));
}

template<typename BaseShape>
std::string OrderedShape2_1<BaseShape>::toString() const {
    return this->underlyingShape->toString();
}

template<typename BaseShape>
const Vector<2> &OrderedShape2_1<BaseShape>::getPosition() const {
    return this->underlyingShape->getPosition();
}

template<typename BaseShape>
void OrderedShape2_1<BaseShape>::translate(const Vector<2> &v) {
    this->underlyingShape->translate(v);
}