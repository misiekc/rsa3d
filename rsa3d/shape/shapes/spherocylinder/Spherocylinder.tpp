/*
 * Spherocylinder.cpp
 *
 *  Created on: Mar 15, 2019
 *      Author: Michal Ciesla
 */

#include "OwnOverlapSC.h"
#include "Stolen2DOverlapSC.h"
#include "../../OrderParameters.h"

static const double EPSILON = 0.0000000001;

template<unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::radius;
template<unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::length;

template <unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::gamma(unsigned short dim){
	double result;
	if(dim % 2==0){
		result = Spherocylinder<DIMENSION>::g20;
		for (unsigned short i=4; i<=dim; i+=2){
			result = (i/2.0)*result;
		}
	}else{ // d%2==1
		result = Spherocylinder<DIMENSION>::g15;
		for (unsigned short i=3; i<=dim; i+=2){
			result = (i/2.0)*result;
		}
	}
	return result;
}

template <unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::volume() {
	double v1 = pow(Spherocylinder<DIMENSION>::radius, DIMENSION)*pow(M_PI, DIMENSION/2.0) / Spherocylinder<DIMENSION>::gamma(DIMENSION);
	double v2 = Spherocylinder<DIMENSION>::length * pow(Spherocylinder<DIMENSION>::radius, DIMENSION-1)*pow(M_PI, (DIMENSION-1)/2.0) / Spherocylinder<DIMENSION>::gamma(DIMENSION-1);
	return v1 + v2;
}


template<unsigned short DIMENSION>
void Spherocylinder<DIMENSION>::initClass(const std::string &attr) {
    std::istringstream attrStream(attr);
    double ratio;

    attrStream >> ratio;
    if (!attrStream)    throw std::runtime_error("Wrong attr format");
    Spherocylinder<DIMENSION>::radius = 1.0;
    Spherocylinder<DIMENSION>::length = 2*Spherocylinder<DIMENSION>::radius*(ratio - 1.0);
    double scaleFactor = pow(Spherocylinder<DIMENSION>::volume(), 1.0/DIMENSION);
    Spherocylinder<DIMENSION>::radius /= scaleFactor;
    Spherocylinder<DIMENSION>::length /= scaleFactor;

    double circumsphareRadius = Spherocylinder<DIMENSION>::radius + Spherocylinder<DIMENSION>::length / 2;
    double insphereRadius = Spherocylinder<DIMENSION>::radius;

    ShapeStaticInfo<DIMENSION, 0> shapeInfo;
    shapeInfo.setCircumsphereRadius(circumsphareRadius);
    shapeInfo.setInsphereRadius(insphereRadius);
    shapeInfo.setExclusionZoneMaxSpan(insphereRadius + circumsphareRadius);
    shapeInfo.setCreateShapeImpl([](RND *rnd) -> Shape<DIMENSION, 0>* {
        if constexpr (DIMENSION == 2)
            return new Spherocylinder<DIMENSION>(Matrix<DIMENSION, DIMENSION>::rotation(
                    2 * M_PI * rnd->nextValue()));
        else if constexpr (DIMENSION == 3)
            return new Spherocylinder<DIMENSION>(Matrix<DIMENSION, DIMENSION>::rotation(
                    2 * M_PI * rnd->nextValue(),
                    std::asin(2 * rnd->nextValue() - 1),
                    2 * M_PI * rnd->nextValue()));
        else
            throw std::runtime_error("Spherocylinder is currently supported only in 2D and 3D.");
    });

    Shape<DIMENSION, 0>::setShapeStaticInfo(shapeInfo);

    // Calculate SpheroCylinder2D params for Stolen2DOverlapSC
	if constexpr (DIMENSION == 2)
        SpheroCylinder2D::calculateStatic(attr);
}

template<unsigned short DIMENSION>
Spherocylinder<DIMENSION>::Spherocylinder(const Matrix<DIMENSION, DIMENSION> &orientation) : orientation(orientation){
}

template<unsigned short DIMENSION>
Vector<DIMENSION> Spherocylinder<DIMENSION>::getEnd(short beginOrEnd) const {
	Vector<DIMENSION> result;
	result[0] = 1.0;
	for (unsigned short i = 1; i < DIMENSION; i++)
		result[i] = 0.0;
	result = (this->getPosition())
			+ (this->orientation * result) * (0.5 * beginOrEnd * Spherocylinder<DIMENSION>::length);
	return result;
}

// Based on
// Copyright 2001 softSurfer, 2012 Dan Sunday
// This code may be freely used, distributed and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
template<unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::distanceFrom(const Spherocylinder *s) const{
	Vector<DIMENSION> t1 = this->getEnd(-1);
	Vector<DIMENSION> t2 = this->getEnd(1);
	Vector<DIMENSION> s1 = s->getEnd(-1);
	Vector<DIMENSION> s2 = s->getEnd(1);

	Vector<DIMENSION> u = t2 - t1;
	Vector<DIMENSION> v = s2 - s1;
	Vector<DIMENSION> w = t1 - s1;
	double a = u * u, b = u * v, c = v * v, d = u * w, e = v * w;

	double D = a * c - b * b;        // always >= 0
	double sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
	double tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

	// compute the line parameters of the two closest points
	if (D < EPSILON) { // the lines are almost parallel
		sN = 0.0;      // force using point P0 on segment S1
		sD = 1.0;      // to prevent possible division by 0.0 later
		tN = e;
		tD = c;
	} else {           // get the closest points on the infinite lines
		sN = (b * e - c * d);
		tN = (a * e - b * d);
		if (sN < 0.0) { // sc < 0 => the s=0 edge is visible
			sN = 0.0;
			tN = e;
			tD = c;
		} else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
			sN = sD;
			tN = e + b;
			tD = c;
		}
	}

	if (tN < 0.0) { // tc < 0 => the t=0 edge is visible
		tN = 0.0;
		// recompute sc for this edge
		if (-d < 0.0)
			sN = 0.0;
		else if (-d > a)
			sN = sD;
		else {
			sN = -d;
			sD = a;
		}
	} else if (tN > tD) { // tc > 1  => the t=1 edge is visible
		tN = tD;
		// recompute sc for this edge
		if ((-d + b) < 0.0)
			sN = 0;
		else if ((-d + b) > a)
			sN = sD;
		else {
			sN = (-d + b);
			sD = a;
		}
	}
	// finally do the division to get sc and tc
	sc = (abs(sN) < EPSILON ? 0.0 : sN / sD);
	tc = (abs(tN) < EPSILON ? 0.0 : tN / tD);

	// get the difference of the two closest points
	Vector<DIMENSION> dP = w + (sc * u) - (tc * v);  // =  S1(sc) - S2(tc)

	return dP.norm();   // return the closest distance
}
//===================================================================

template<unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::distanceFrom(Vector<DIMENSION> *v) const{
	Vector<DIMENSION> t1 = this->getEnd(-1);
	Vector<DIMENSION> t2 = this->getEnd(1);

	Vector<DIMENSION> a = t2 - t1;
	Vector<DIMENSION> b = *v - t1;
	Vector<DIMENSION> c = *v - t2;

	if ( ((a*b) * (a*c))>0 ){
		return std::min(b.norm(), c.norm());
	}else{
		return (b - a*(a*b)).norm();
	}
}


template<unsigned short DIMENSION>
bool Spherocylinder<DIMENSION>::overlap(BoundaryConditions<DIMENSION> *bc, const Shape<DIMENSION, 0> *s) const {
    using EarlyRejectionResult = typename Shape<DIMENSION, 0>::EarlyRejectionResult;
    switch (this->overlapEarlyRejection(bc, s)) {
        case EarlyRejectionResult::TRUE:      return true;
        case EarlyRejectionResult::FALSE:     return false;
        case EarlyRejectionResult::UNKNOWN:   break;
    }

    Spherocylinder<DIMENSION> other = dynamic_cast<const Spherocylinder<DIMENSION>&>(*s);
    this->applyBC(bc, &other);
	return this->distanceFrom(&other) < 2*Spherocylinder<DIMENSION>::getRadius();
}

template<unsigned short DIMENSION>
bool Spherocylinder<DIMENSION>::pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &position, const Orientation<0> &orientation, double orientationRange) const {
    Vector<DIMENSION> bcPos = position + bc->getTranslation(this->getPosition(), position);

    double d = this->distanceFrom(&bcPos);
    return (d <= 2.0*Spherocylinder<DIMENSION>::radius);
}

template<unsigned short DIMENSION>
std::vector<double> Spherocylinder<DIMENSION>::calculateOrder(const OrderCalculable *other) const {
    auto &s = dynamic_cast<const Spherocylinder<DIMENSION>&>(*other);    // use reference so dynamic_cast would throw on error

	Vector<DIMENSION> thisAxis = this->getEnd(1) - this->getEnd(-1);
	Vector<DIMENSION> otherAxis = s.getEnd(1) - s.getEnd(-1);

    return { legendre::P2(thisAxis*otherAxis) };
}


template<unsigned short DIMENSION>
void Spherocylinder<DIMENSION>::store(std::ostream &f) const {
    Shape<DIMENSION, 0>::store(f);
    double d;
    for (unsigned short i=0; i<DIMENSION; i++){
        for (unsigned short j=0; j<DIMENSION; j++){
            d = this->orientation(i, j);
            f.write((char *)(&d), sizeof(double));
        }
    }
}

template<unsigned short DIMENSION>
void Spherocylinder<DIMENSION>::restore(std::istream &f) {
    Shape<DIMENSION, 0>::restore(f);
    double d;
    for (unsigned short i=0; i<DIMENSION; i++){
        for (unsigned short j=0; j<DIMENSION; j++){
            f.read((char *)&d, sizeof(double));
            this->orientation(i, j) = d;
        }
    }
}

template<unsigned short DIMENSION>
Shape<DIMENSION, 0> *Spherocylinder<DIMENSION>::clone() const {
    return new Spherocylinder(*this);
}

template<unsigned short DIMENSION>
std::string Spherocylinder<DIMENSION>::toPovray() const {
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);
	Vector<DIMENSION> t1 = this->getEnd(-1);
	Vector<DIMENSION> t2 = this->getEnd(1);

	if constexpr (DIMENSION==2) {
        out << "  cylinder { < " << t1[0] << ", " << t1[1] << ", 0.0 >, < " << t2[0] << ", " << t2[1] << ", 0.0 >, " << Spherocylinder<DIMENSION>::radius << std::endl;
        out << "             texture { pigment { color Red } }" << std::endl;
        out << "  }" << std::endl;
        out << "  sphere { < " << t1[0] << ", " << t1[1] << ", 0.0 >, " << Spherocylinder<DIMENSION>::radius << std::endl;
        out << "             texture { pigment { color Red } }" << std::endl;
        out << "  }" << std::endl;
        out << "  sphere { < " << t2[0] << ", " << t2[1] << ", 0.0 >, " << Spherocylinder<DIMENSION>::radius << std::endl;
        out << "             texture { pigment { color Red } }" << std::endl;
        out << "  }" << std::endl;
    } else if constexpr (DIMENSION==3) {
        out << "  cylinder { < " << t1[0] << ", " << t1[1] << ", " << t1[2] << ">, < " << t2[0] << ", " << t2[1] << ", " << t2[2] << ">, " << Spherocylinder<DIMENSION>::radius << std::endl;
        out << "             texture { pigment { color Red } }" << std::endl;
        out << "  }" << std::endl;
        out << "  sphere { < " << t1[0] << ", " << t1[1] << ", " << t1[2] << ">, " << Spherocylinder<DIMENSION>::radius << std::endl;
        out << "             texture { pigment { color Red } }" << std::endl;
        out << "  }" << std::endl;
        out << "  sphere { < " << t2[0] << ", " << t2[1] << ", " << t2[2] << ">, " << Spherocylinder<DIMENSION>::radius << std::endl;
        out << "             texture { pigment { color Red } }" << std::endl;
        out << "  }" << std::endl;
    } else {
        throw std::runtime_error("Spherocylinder::toPovray() supported only for 2D and 3D");
    }
	return out.str();
}

template<unsigned short DIMENSION>
std::string Spherocylinder<DIMENSION>::toWolfram() const {
    std::stringstream out;
    out << std::fixed;

	if constexpr (DIMENSION == 2) {
        out << "GeometricTransformation[{Rectangle[{-" << length / 2 << ", -" << radius << "}, {" << length / 2 << ", "
            << radius << "}]," << std::endl;
        out << "    Disk[{-" << length / 2 << ", 0}, " << radius << "]," << std::endl;
        out << "    Disk[{" << length / 2 << ", 0}, " << radius << "]}," << std::endl;
        out << "    {{{" << this->orientation(0, 0) << ", " << this->orientation(0, 1) << "}, ";
        out << "{" << this->orientation(1, 0) << ", " << this->orientation(1, 1) << "}}, ";
        out << this->getPosition() << "}]";
    } else if (DIMENSION == 3) {
        Vector<3> beg = this->getEnd(-1);
        Vector<3> end = this->getEnd(1);
        out << "CapsuleShape[{" << beg << ", " << end << "}, " << Spherocylinder::radius << "]";
    } else {
        throw std::runtime_error("Spherocylinder::toWolfram() supported only for 2D and 3D");
    }

    return out.str();
}

template<unsigned short DIMENSION>
std::vector<std::string> Spherocylinder<DIMENSION>::getSupportedStrategies() const {
    if constexpr (DIMENSION == 2)
        return {"own", "stolen_from_sc2d"};
    else
        return {"own"};
}

template<unsigned short DIMENSION>
OverlapStrategy<DIMENSION, 0> *Spherocylinder<DIMENSION>::createStrategy(const std::string &name) const {
	if constexpr (DIMENSION == 2)
        if (name == "stolen_from_sc2d")
            return new Stolen2DOverlapSC;

    if (name == "own")
        return new OwnOverlapSC<DIMENSION>;

    return nullptr;
}

template<unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::getLength() {
    return length;
}

template<unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::getRadius() {
    return radius;
}

template<unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::getAngle() const {
    return std::atan2(this->orientation(1, 0), this->orientation(0, 0));
}
