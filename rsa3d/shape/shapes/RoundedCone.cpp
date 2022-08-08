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

double RoundedCone::R;
double RoundedCone::r;
double RoundedCone::l;

RoundedCone::RoundedCone(const Matrix<3, 3> &orientation) : GenericXCShape(orientation)
{
}

RoundedCone::~RoundedCone() {
}

// needed to be chcecked
double RoundedCone::volume(){
	return
			(M_PI*
					(l*l*l*l*(r*r + r*R + R*R) +
					2*l*l*(r - R)*(r - R)*(r*r + r*R + R*R) -
					(r - R)*(r - R)*(r - R)*(r - R)*(r*r + r*R + R*R) +
					2*l*l*l*(r*r*r + R*R*R))
			)/(3*l*l*l);
}

/**
 * @param attr two parameters:
 * - radius of the smaller cap
 * - distance between centers of spheres
 */
void RoundedCone::initClass(const std::string &attr){
    std::istringstream attrStream(attr);
    RoundedCone::R=1.0;
    attrStream >> RoundedCone::r >> RoundedCone::l;
    if (!attrStream)    throw ValidationException("Wrong attr format");
    ValidateMsg(r > 0 && l > 0, "RoundedCone parameters should be positive");
    ValidateMsg(r <= 1, "RoundedCone smaller sphere radius should be smaller than 1.0");

    double v = volume();
    double factor = std::pow(v, 1.0/3.0);
    R /= factor;
    r /= factor;
    l /= factor;
    v = volume();

    BodyBuilder bb;
    bb.sphere(R);
    bb.move(0, 0, -l/2);
    bb.sphere(r);
    bb.move(0, 0, l/2);
    bb.wrap();
    GenericXCShape::shapeModel = bb.getCollideGeometry();

    bb.clear();
    double d = l*(R+r)/std::sqrt(4*l*l + (R-r)*(R-r));

    bb.sphere(R+d);
    bb.move(0, 0, -l/2);
    bb.sphere(r+d);
    bb.move(0, 0, l/2);
    bb.wrap();
    GenericXCShape::evModel = bb.getCollideGeometry();


    ShapeStaticInfo<3, 0> shapeInfo;
    shapeInfo.setCircumsphereRadius(0.5*l+R);
    shapeInfo.setInsphereRadius(d);
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

// not supported yet
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


