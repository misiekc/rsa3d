/*
 * GenericXCShape.cpp
 *
 *  Created on: Aug 7, 2022
 *      Author: ciesla
 */

#include "GenericXCShape.h"
#include "../../geometry/xenocollide/MapPtr.h"
#include "../../geometry/xenocollide/BodyBuilder.h"
#include "../../geometry/xenocollide/Collide.h"

MapPtr<CollideGeometry> GenericXCShape::shapeModel;
MapPtr<CollideGeometry> GenericXCShape::evModel;

GenericXCShape::GenericXCShape(const Matrix<3, 3> &orientation) : orientation(orientation) {
}

GenericXCShape::~GenericXCShape() {
}

/**
 * @param attr parameters are:
 * - radius of the insphere
 * - radius of the circumsphere
 * - commands that builds shape and its excluded volume.
 * Commands are interpreted by BodyBuilder class and should be divided by '&' character.
 * Program expects two shapes. If there is only one shape defined it also takes role of its excluded volume
 * Example of definition of a single shape
 * 0.5 1.0 sphere 0.5 & move 0 0 -0.5 & sphere 1 & move 0 0 0.5 & wrap
 */
void GenericXCShape::initClass(const std::string &attr){
	std::stringstream ss(attr);
	double R, r;
	ss >> r >> R;
	std::string commands;
	std::getline(ss, commands, '\0');
	std::string script = "script " + commands;

	BodyBuilder bb;
    bb.ProcessCommand(script);
	GenericXCShape::evModel = bb.getCollideGeometry();
    if (bb.getModelStackSize()==2){
    	bb.pop();
    }
	GenericXCShape::shapeModel = bb.getCollideGeometry();

    ShapeStaticInfo<3, 0> shapeInfo;
    shapeInfo.setCircumsphereRadius(R);
    shapeInfo.setInsphereRadius(r);
    shapeInfo.setCreateShapeImpl([](RND *rnd) -> Shape* {
        return new GenericXCShape(Matrix<3, 3>::rotation(
                2 * M_PI * rnd->nextValue(),
                std::asin(2 * rnd->nextValue() - 1),
                2 * M_PI * rnd->nextValue()));
    });
    Shape::setShapeStaticInfo(shapeInfo);
}

Shape<3, 0> *GenericXCShape::clone() const {
    return new GenericXCShape(*this);
}

bool GenericXCShape::overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const {
    switch (this->overlapEarlyRejection(bc, s)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }

	GenericXCShape cone = dynamic_cast<const GenericXCShape &>(*s);
    this->applyBC(bc, &cone);
    Quat q1(this->orientation);
    Quat q2(cone.orientation);
    bool result = Collide::Intersect(*(shapeModel), q1, this->getPosition(), *(shapeModel), q2, cone.getPosition(), 1.0e-12);
	return result;
}

bool GenericXCShape::pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation, double orientationRange) const {
    Vector<3> bcPos = position + bc->getTranslation(this->getPosition(), position);
    CollideGeometry *point = new CollidePoint(Vector<3>({0, 0, 0}));
    Quat q1(this->orientation);
    Quat q2(Matrix<3,3>::identity());
    bool result = Collide::Intersect(*(evModel), q1, this->getPosition(), *point, q2, bcPos, 1.0e-12);
    return result;
}

void GenericXCShape::store(std::ostream &f) const {
    Shape::store(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            d = this->orientation(i, j);
            f.write((char *)(&d), sizeof(double));
        }
    }
}

void GenericXCShape::restore(std::istream &f) {
    Shape::restore(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            f.read((char *)&d, sizeof(double));
            this->orientation(i, j) = d;
        }
    }
}


