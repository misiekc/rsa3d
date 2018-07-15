/*
 * Sphere.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#include <math.h>
#include <iostream>

#include "../../BoundaryConditions.h"


template <unsigned short DIMENSION>
double Sphere<DIMENSION>::radius;

template <unsigned short DIMENSION>
Sphere<DIMENSION>::~Sphere() {
}

template <unsigned short DIMENSION>
double Sphere<DIMENSION>::gamma(){
	double result;
	if(DIMENSION % 2==0){
		result = Sphere::g20;
		for (unsigned short i=4; i<=DIMENSION; i+=2){
			result = (i/2.0)*result;
		}
	}else{ // d%2==1
		result = Sphere::g15;
		for (unsigned short i=3; i<=DIMENSION; i+=2){
			result = (i/2.0)*result;
		}
	}
	return result;
}

template <unsigned short DIMENSION>
double Sphere<DIMENSION>::volume() {
	return pow(M_PI, DIMENSION/2.0) / gamma();
}

template <unsigned short DIMENSION>
void Sphere<DIMENSION>::initClass(const std::string &args){
	Sphere<DIMENSION>::radius = pow(1.0/Sphere::volume(), 1.0/DIMENSION);
	Shape<DIMENSION, 0>::setNeighbourListCellSize(2.0*Sphere::radius);
	Shape<DIMENSION, 0>::setVoxelSpatialSize(pow(2.0, 1.0/DIMENSION)*Sphere::radius);
//	Shape<DIMENSION, 0>::setCreateShapeImpl([](RND *rnd) -> Shape<DIMENSION, 0> {
//        return new Sphere<DIMENSION, 0>;
//    });
    Shape<DIMENSION, 0>::template setDefaultCreateShapeImpl <Sphere<DIMENSION>> ();
}

template <unsigned short DIMENSION>
Sphere<DIMENSION>::Sphere() : ConvexShape<DIMENSION, 0>(){
	this->r = Sphere<DIMENSION>::radius;
}

template <unsigned short DIMENSION>
bool Sphere<DIMENSION>::overlap(BoundaryConditions<DIMENSION> *bc, const Shape<DIMENSION, 0> *s) const {
	const Sphere *sd = dynamic_cast<const Sphere<DIMENSION>*>(s);
	double d2 = bc->distance2(this->getPosition(), sd->getPosition());
	double r2 = this->r + sd->r;
	r2 *= r2;

	return (d2 < r2);
}

template <unsigned short DIMENSION>
double Sphere<DIMENSION>::getVolume() const {
	return Sphere<DIMENSION>::volume()*pow(this->r, DIMENSION);
}

template <unsigned short DIMENSION>
bool Sphere<DIMENSION>::pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &da) const {
	double d2;
	if (bc!=NULL)
		d2 = bc->distance2(da, this->getPosition());
	else{
		d2 = 0.0;
		for(unsigned short i=0; i<DIMENSION; i++)
			d2 += (da[i]-this->getPosition()[i])*(da[i]-this->getPosition()[i]);
	}
	return (d2<4.0*this->r*this->r);
}

template <unsigned short DIMENSION>
bool Sphere<DIMENSION>::pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &position,
									const Orientation<0> &orientation, double orientationRange) const {
    return this->pointInside(bc, position);
}

template <unsigned short DIMENSION>
double Sphere<DIMENSION>::minDistance(Shape<DIMENSION, 0> *s) const{
	return 2.0*this->radius;
}

template <unsigned short DIMENSION>
std::string Sphere<DIMENSION>::toPovray() const{
	std::string s;
    Vector<DIMENSION> position = this->getPosition();
	double r0 = this->r - 0.01;

	if (DIMENSION==2){
		s = "  disc { < ";
		for(unsigned short i=0; i<DIMENSION; i++)
			s += std::to_string(position[i]) + ", ";
		s += "0.05>, <0.0, 0.0, 1.0>, " + std::to_string(this->r) +"\r\n    texture { pigment { color Red } }\r\n  }\r\n";
/*
		s += "  disc { < ";
		for(unsigned short i=0; i<DIMENSION; i++)
			s += std::to_string(this->position[i]) + ", ";
		s += "0.01>, <0.0, 0.0, 1.0>, " + std::to_string(2*r0) + "\r\n    texture { pigment { color Yellow } }\r\n  }\r\n";
*/
		s += "  disc { < ";
		for(unsigned short i=0; i<DIMENSION; i++)
			s += std::to_string(position[i]) + ", ";
		s += "0.03>, <0.0, 0.0, 1.0>, " + std::to_string(2*this->r) + ", " + std::to_string(2*r0) + "\r\n    texture { pigment { color Black } }\r\n  }\r\n";

/*
		s += "  text { ttf \"timrom.ttf\" \"" + std::to_string(this->no) + "\" 1, 0 pigment { color Black } scale 0.5 translate < ";
		for(unsigned char i=0; i<this->dimension; i++)
			s += std::to_string(this->position[i]) + ", ";
		s +=  "0.0003> }\r\n";
*/
	}else{
		s = "  sphere { < ";
		for(unsigned short i=0; i<DIMENSION; i++)
			s += std::to_string(position[i]) + ", ";
		s += "0.0>, " + std::to_string(this->r) +"\n    texture { pigment { color Red } }\r\n  }\r\n";
	}
	return s;
}

template <unsigned short DIMENSION>
void Sphere<DIMENSION>::store(std::ostream &f) const{
	Shape<DIMENSION, 0>::store(f);
	f.write((char *)(&this->r), sizeof(double));
}

template <unsigned short DIMENSION>
void Sphere<DIMENSION>::restore(std::istream &f){
	Shape<DIMENSION, 0>::restore(f);
	f.read((char *)(&this->r), sizeof(double));
}

template<unsigned short DIMENSION>
Shape<DIMENSION, 0> *Sphere<DIMENSION>::clone() const {
    return new Sphere(*this);
}
