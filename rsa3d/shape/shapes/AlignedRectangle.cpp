//
// Created by Michal Ciesla on 8/06/2022.
//

#include <vector>
#include "AlignedRectangle.h"
#include <cmath>

std::vector<double> AlignedRectangle::allowedOrientations;

AlignedRectangle::AlignedRectangle() : Rectangle() {
}

void AlignedRectangle::initClass(const std::string &args) {

	Rectangle::initClass(args);
	std::istringstream in(args);
	double ratio, orientationsNo;
	in >> ratio;
	in >> orientationsNo;
	double da = M_PI / orientationsNo;
	for (double angle=0; angle<M_PI; angle += da){
		AlignedRectangle::allowedOrientations.push_back(angle);
	}

    ShapeStaticInfo<2, 1> shapeInfo = Shape::getShapeStaticInfo();
    shapeInfo.setDefaultCreateShapeImpl<AlignedRectangle>();

    Shape::setShapeStaticInfo(shapeInfo);
}

void AlignedRectangle::setAngle(double angle) {
	if (angle<0) angle += M_PI;
	if (angle>M_PI) angle -= M_PI;
	size_t minIndex = -1;
	double dist, min = std::numeric_limits<double>::infinity();

	for(size_t index=0; index<AlignedRectangle::allowedOrientations.size(); index++){
		dist = std::min(
			std::abs(angle - AlignedRectangle::allowedOrientations[index]),
			std::abs(angle - AlignedRectangle::allowedOrientations[index]-M_PI)
			);
		if (dist<min){
			min = dist;
			minIndex = index;
		}
	}
	Rectangle::setAngle(AlignedRectangle::allowedOrientations[minIndex]);
}


bool AlignedRectangle::pointInside(BoundaryConditions<2> *bc, const Vector<2> &da, double angleFrom, double angleTo) const {
	// let's have it ordered
    if (angleFrom > angleTo) {
        double tmp = angleFrom;
        angleFrom = angleTo;
        angleTo = tmp;
    }

    bool allowedAnglesInRange = false;
	for(double angle: AlignedRectangle::allowedOrientations){
		if (angle>=angleFrom && angle<=angleTo){
			allowedAnglesInRange = true;
			break;
		}
	}

	if (!allowedAnglesInRange)
		return true;
	return Rectangle::pointInside(bc, da, angleFrom, angleTo);
}

