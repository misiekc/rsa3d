/*
 * Spherocylinder.cpp
 *
 *  Created on: Mar 15, 2019
 *      Author: Michal Ciesla
 */

#include "../../Vector.h"

static const double EPSILON = 0.0000000001;

template<unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::radius;
template<unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::lenght;


template<unsigned short DIMENSION>
void Spherocylinder<DIMENSION>::initClass(const std::string &attr) {
    std::istringstream attrStream(attr);
    double ratio;

    attrStream >> ratio;
    if (!attrStream)    throw std::runtime_error("Wrong attr format");
    Spherocylinder<DIMENSION>::radius = std::pow(M_PI*(2.0*(ratio - 1.0) + 4.0/3.0), -1.0/3.0);
    Spherocylinder<DIMENSION>::lenght = 2*Spherocylinder<DIMENSION>::radius*(ratio - 1.0);
    Shape<DIMENSION, 0>::setNeighbourListCellSize(2.0*Spherocylinder<DIMENSION>::radius + Spherocylinder<DIMENSION>::lenght);
    Shape<DIMENSION, 0>::setVoxelSpatialSize(M_SQRT2 * Spherocylinder<DIMENSION>::radius);

    Shape<DIMENSION, 0>::setCreateShapeImpl([](RND *rnd) -> Shape<DIMENSION, 0>* {
        return new Spherocylinder<DIMENSION>(Matrix<DIMENSION, DIMENSION>::rotation(
                2 * M_PI * rnd->nextValue(),
                std::asin(2 * rnd->nextValue() - 1),
                2 * M_PI * rnd->nextValue()));
    });
}

template<unsigned short DIMENSION>
Spherocylinder<DIMENSION>::Spherocylinder(const Matrix<DIMENSION, DIMENSION> &orientation) : orientation(orientation){
}

template<unsigned short DIMENSION>
void Spherocylinder<DIMENSION>::getEnd(short beginOrEnd, Vector<DIMENSION> *result) const {
	(*result)[0] = 1.0;
	for (unsigned short i = 1; i < DIMENSION; i++)
		(*result)[i] = 0.0;
	*result = (this->getPosition())
			+ (this->orientation * (*result)) * (0.5 * beginOrEnd * Spherocylinder<DIMENSION>::lenght);
}

// Based on
// Copyright 2001 softSurfer, 2012 Dan Sunday
// This code may be freely used, distributed and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
template<unsigned short DIMENSION>
double Spherocylinder<DIMENSION>::distanceFrom(Spherocylinder *s) const{
	Vector<DIMENSION> t1, t2, s1, s2;
	this->getEnd(-1, &t1);
	this->getEnd(1, &t2);
	s->getEnd(-1, &s1);
	s->getEnd(1, &s2);

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
	Vector<DIMENSION> t1, t2;
	this->getEnd(-1, &t1);
	this->getEnd(1, &t2);

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
    Spherocylinder<DIMENSION> other = dynamic_cast<const Spherocylinder<DIMENSION>&>(*s);
    this->applyBC(bc, &other);
    return this->distanceFrom(&other) < 2*Spherocylinder<DIMENSION>::radius;
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

	Vector<DIMENSION> t1, t2, s1, s2;
	this->getEnd(-1, &t1);
	this->getEnd(1, &t2);
	s.getEnd(-1, &s1);
	s.getEnd(1, &s2);

    std::vector<double> result(1);
    result[0] = P2((t2-t1)*(s2-s1));
    return result;
}


template<unsigned short DIMENSION>
void Spherocylinder<DIMENSION>::store(std::ostream &f) const {
    Shape<DIMENSION, 0>::store(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            d = this->orientation(i, j);
            f.write((char *)(&d), sizeof(double));
        }
    }
}

template<unsigned short DIMENSION>
void Spherocylinder<DIMENSION>::restore(std::istream &f) {
    Shape<DIMENSION, 0>::restore(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
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
	Vector<DIMENSION> t1, t2;
	this->getEnd(-1, &t1);
	this->getEnd(1, &t2);
	out << "  cylinder { < " << t1[0] << ", " << t1[1] << ", " << t1[2] << ">, < " << t2[0] << ", " << t2[1] << ", " << t2[2] << ">, " << Spherocylinder<DIMENSION>::radius << std::endl;
	out << "             texture { pigment { color Red } }" << std::endl;
	out << "  }" << std::endl;
	out << "  sphere { < " << t1[0] << ", " << t1[1] << ", " << t1[2] << ">, " << Spherocylinder<DIMENSION>::radius << std::endl;
	out << "             texture { pigment { color Red } }" << std::endl;
	out << "  }" << std::endl;
	out << "  sphere { < " << t2[0] << ", " << t2[1] << ", " << t2[2] << ">, " << Spherocylinder<DIMENSION>::radius << std::endl;
	out << "             texture { pigment { color Red } }" << std::endl;
	out << "  }" << std::endl;
	return out.str();
}
