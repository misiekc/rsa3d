/*
 * RoundedCone.cpp
 *
 *  Created on: Aug 5, 2022
 *      Author: ciesla
 */

#include "RoundedCone.h"
#include "../../geometry/xenocollide/BodyBuilder.h"
#include "../../geometry/xenocollide/Collide.h"
#include "../../geometry/xenocollide/Quat.h"

MapPtr<CollideGeometry> RoundedCone::shapeModel;
MapPtr<CollideGeometry> RoundedCone::evModel;
double RoundedCone::R;
double RoundedCone::r;
double RoundedCone::l;

RoundedCone::RoundedCone(const Matrix<3, 3> &orientation) : orientation(orientation)
{

}


RoundedCone::~RoundedCone() {
	// TODO Auto-generated destructor stub
}

void RoundedCone::initClass(const std::string &attr){
    std::istringstream attrStream(attr);
    RoundedCone::R=1.0;
    attrStream >> RoundedCone::r >> RoundedCone::l;
    if (!attrStream)    throw ValidationException("Wrong attr format");
    ValidateMsg(r > 0 && l > 0, "RoundedCone parameters should be positive");
    ValidateMsg(r <= 1, "RoundedCone smaller sphere radius should be smaller than 1.0");

    double v = (1.0/3.0) * M_PI *
    		(2*R*R*R + 2*r*r*r + (l/(R-r))*(R*R*R-r*r*r));
    double factor = std::pow(v, 1.0/3.0);
    R /= factor;
    r /= factor;
    l /= factor;

    BodyBuilder bb;
    bb.sphere(R);
    bb.move(0, 0, -l/2);
    bb.sphere(r);
    bb.move(0, 0, l/2);
    bb.wrap();
    RoundedCone::shapeModel = bb.getCollideGeometry();

    bb.clear();

    bb.sphere(R+r);
    bb.move(0, 0, -l/2);
    bb.sphere(r+r);
    bb.move(0, 0, l/2);
    bb.wrap();
    RoundedCone::evModel = bb.getCollideGeometry();


    ShapeStaticInfo<3, 0> shapeInfo;
    shapeInfo.setCircumsphereRadius(l+R+r);
    shapeInfo.setInsphereRadius(r);
//    shapeInfo.setExclusionZoneMaxSpan(a + c);
    shapeInfo.setCreateShapeImpl([](RND *rnd) -> Shape* {
        return new RoundedCone(Matrix<3, 3>::rotation(
                2 * M_PI * rnd->nextValue(),
                std::asin(2 * rnd->nextValue() - 1),
                2 * M_PI * rnd->nextValue()));
    });
    Shape::setShapeStaticInfo(shapeInfo);
}

Shape<3, 0> *RoundedCone::clone() const {
    return new RoundedCone(*this);
}

bool RoundedCone::overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const {
//    switch (this->overlapEarlyRejection(bc, s)) {
//        case TRUE:      return true;
//        case FALSE:     return false;
//        case UNKNOWN:   break;
//    }

    RoundedCone cone = dynamic_cast<const RoundedCone &>(*s);
    this->applyBC(bc, &cone);
    Quat q1(this->orientation);
    Quat q2(cone.orientation);
//    CollideSphere geom(1.0);
//    geom = *model;
//	bool result = Intersect(geom, q1, this->getPosition(), geom, q2, cone.getPosition(), 0.1);
//	double dist = (this->getPosition()-cone.getPosition()).norm();
    bool result = Collide::Intersect(*(shapeModel), q1, this->getPosition(), *(shapeModel), q2, cone.getPosition(), 1.0e-12);
//	if (result==false && dist<1.0){
//		std::cout << this->getPosition() << ", " << cone.getPosition() << std::endl;
//	}
	return result;
}

bool RoundedCone::pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation, double orientationRange) const {
    // Transform point coordinates to Ellipsoid coordinate system
    Vector<3> bcPos = position + bc->getTranslation(this->getPosition(), position);
    CollideGeometry *point = new CollidePoint(Vector<3>({0, 0, 0}));
    Quat q1(this->orientation);
    Quat q2(Matrix<3,3>::identity());
    bool result = Collide::Intersect(*(evModel), q1, this->getPosition(), *point, q2, bcPos, 1.0e-12);
    return result;
}

void RoundedCone::store(std::ostream &f) const {
    Shape::store(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            d = this->orientation(i, j);
            f.write((char *)(&d), sizeof(double));
        }
    }
}

void RoundedCone::restore(std::istream &f) {
    Shape::restore(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            f.read((char *)&d, sizeof(double));
            this->orientation(i, j) = d;
        }
    }
}

std::string RoundedCone::toPovray() const {
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);
    Vector<3> position = this->getPosition();
	out << "  sphere { < 0, 0, 0 >, 1.0" << std::endl;
	out << "    matrix <" << std::endl;
	out << "      " << this->orientation(0, 0) << ", " << this->orientation(1, 0) << ", " << this->orientation(2, 0) << ", " << std::endl;
    out << "      " << this->orientation(0, 1) << ", " << this->orientation(1, 1) << ", " << this->orientation(2, 1) << ", " << std::endl;
    out << "      " << this->orientation(0, 2) << ", " << this->orientation(1, 2) << ", " << this->orientation(2, 2) << ", " << std::endl;
    out << "      " << position[0] << ", " << position[1] << ", " << position[2] << std::endl;
    out << "    >" << std::endl;
	out << "    texture { pigment { color Red } }" << std::endl;
	out << "  }" << std::endl;
	return out.str();
}

std::string RoundedCone::toWolfram() const {
	int ballsNo = 10;
    std::stringstream out;
    out << std::fixed;
    out << "GeometricTransformation[" << std::endl;
    out << "    {";
    for(int i=0; i<ballsNo; i++){
    	out << "Sphere[{0, 0, "<< -l/2 + l*(((double)i)/(ballsNo-1)) << "}, " << ((ballsNo-1.0-i)*R + i*r)/(ballsNo-1) <<"]";
    	if (i<ballsNo-1)
    		out <<", ";
    }
    out << "}, " << std::endl;
    out << "    AffineTransform[" << std::endl;
    out << "        {{{" << this->orientation(0, 0) << ", " << this->orientation(0, 1) << ", " << this->orientation(0, 2) << "}," << std::endl;
    out << "          {" << this->orientation(1, 0) << ", " << this->orientation(1, 1) << ", " << this->orientation(1, 2) << "}," << std::endl;
    out << "          {" << this->orientation(2, 0) << ", " << this->orientation(2, 1) << ", " << this->orientation(2, 2) << "}}," << std::endl;
    out << "          " << this->getPosition() << "}]]";
    return out.str();
}


