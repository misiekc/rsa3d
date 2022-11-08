//
// Created by Michal Ciesla on 20.10.2022.
//

#include <sstream>
#include <vector>

#include "DiscreteOrientationsShape2_1.h"

std::vector<RSAOrientation> DiscreteOrientationsShape2_1::allowedOrientations{};

/* Arguments in format: (angle distribution average) (angle distribution sigma) (encapsulated class parameters) */
void DiscreteOrientationsShape2_1::initClass(const std::string &args, InitClassFunction baseShapeInitClassFunction) {
    std::istringstream argsStream(args);
    size_t n;
    argsStream >> n;
    ValidateMsg(n>0, "At least one angle needed");
    std::string sToken;
    argsStream >> sToken;
    if (sToken=="auto"){
    	for(double angle = 0; angle<2.0*M_PI; angle +=2.0*M_PI/n){
        	RSAOrientation orientation({angle});
        	allowedOrientations.push_back(orientation);
    	}
    }else{
    	for(size_t i=0; i<n; i++){
    		double decAngle = std::atoi(sToken.c_str());
    		ValidateMsg(decAngle >= 0 && decAngle < 360, "Angle should be from [0, 360) interval");
    		RSAOrientation orientation({decAngle*M_PI/180.0});
    		allowedOrientations.push_back(orientation);
    		if (i<n-1){
    			argsStream >> sToken;
    		    ValidateMsg(argsStream, " n angle1 angle 2 ... angleN in degrees \n or \n n auto");
    		}
    	}
    }

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

Shape<2, 0> *DiscreteOrientationsShape2_1::createShape(RND *rnd) {
    auto createShape = Shape<2, 1>::getCreateShapeImpl();
    auto underlyingShape = std::unique_ptr<Shape<2, 1>>(createShape(rnd));

    Orientation<1> orientation = allowedOrientations[(size_t)(allowedOrientations.size()*rnd->nextValue())];
    underlyingShape->rotate(orientation);

    return new DiscreteOrientationsShape2_1(std::move(underlyingShape));
}

DiscreteOrientationsShape2_1::DiscreteOrientationsShape2_1(std::unique_ptr<Shape<2, 1>> underlyingShape) {
    this->underlyingShape = std::move(underlyingShape);
}

DiscreteOrientationsShape2_1::DiscreteOrientationsShape2_1(const DiscreteOrientationsShape2_1 &other) : Shape(other) {
    this->underlyingShape.reset(other.underlyingShape->clone());
}

DiscreteOrientationsShape2_1 &DiscreteOrientationsShape2_1::operator=(const DiscreteOrientationsShape2_1 &other) {
    Shape::operator=(other);
    this->underlyingShape.reset(other.underlyingShape->clone());
    return *this;
}

bool DiscreteOrientationsShape2_1::overlap(BoundaryConditions<2> *bc, const Shape<2, 0> *s) const {
    auto &orderedOther = dynamic_cast<const DiscreteOrientationsShape2_1&>(*s);
    return this->underlyingShape->overlap(bc, orderedOther.underlyingShape.get());
}

double DiscreteOrientationsShape2_1::getVolume(unsigned short dim) const {
    return this->underlyingShape->getVolume(dim);
}

bool DiscreteOrientationsShape2_1::voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition,
                                  const Orientation<0> &voxelOrientation, double spatialSize,
                                  double angularSize) const {
	for(Orientation<1> orientation: allowedOrientations)
		if (!this->underlyingShape->voxelInside(bc, voxelPosition, orientation, spatialSize, 0.0))
			return false;
	return true;
}

double DiscreteOrientationsShape2_1::minDistance(const Shape<2, 0> *s) const {
    auto &orderedOther = dynamic_cast<const DiscreteOrientationsShape2_1&>(*s);
    return this->underlyingShape->minDistance(orderedOther.underlyingShape.get());
}

std::string DiscreteOrientationsShape2_1::toPovray() const {
    return this->underlyingShape->toPovray();
}

std::string DiscreteOrientationsShape2_1::toWolfram() const {
    return this->underlyingShape->toWolfram();
}

void DiscreteOrientationsShape2_1::store(std::ostream &f) const {
    return this->underlyingShape->store(f);
}

void DiscreteOrientationsShape2_1::restore(std::istream &f) {
    return this->underlyingShape->restore(f);
}

Shape<2, 0> *DiscreteOrientationsShape2_1::clone() const {
    return new DiscreteOrientationsShape2_1(*this);
}

std::string DiscreteOrientationsShape2_1::toString() const {
    return this->underlyingShape->toString();
}

const Vector<2> &DiscreteOrientationsShape2_1::getPosition() const {
    return this->underlyingShape->getPosition();
}

void DiscreteOrientationsShape2_1::translate(const Vector<2> &v) {
    this->underlyingShape->translate(v);
}
