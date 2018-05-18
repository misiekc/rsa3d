/*
 * Voxel.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#include "Voxel.h"

#include <sstream>
#include <stdexcept>
#include <array>
#include <iostream>


Voxel::Voxel(){
	this->lastAnalyzed = 0;
	this->depth = 0;
}


Voxel::Voxel(double* pos, const std::array<double, RSA_ANGULAR_DIMENSION> &angle){
	std::copy(pos, pos+RSA_SPATIAL_DIMENSION, this->position);
	this->orientation = angle;
	this->lastAnalyzed = 0;
	this->depth = 0;
}

bool Voxel::isInside(double *pos, double size){
	for(int i=0; i<RSA_SPATIAL_DIMENSION; i++){
		if (pos[i]<this->position[i])
			return false;
		if (pos[i]>=(this->position[i]+size))
			return false;
	}
	return true;
}


bool Voxel::isInside(double *pos, double size, const std::array<double, RSA_ANGULAR_DIMENSION> &angle, double asize){
	for(int i=0; i<RSA_SPATIAL_DIMENSION; i++){
		if (pos[i]<this->position[i])
			return false;
		if (pos[i]>=(this->position[i]+size))
			return false;
	}

	for(int i=0; i<RSA_ANGULAR_DIMENSION; i++){
		if (angle[i]<this->orientation[i])
			return false;
		if (angle[i]>=(this->orientation[i]+asize))
			return false;
	}

	return true;
}



double* Voxel::getPosition(){
	return this->position;
}


std::array<double, RSA_ANGULAR_DIMENSION> Voxel::getOrientation(){
	return this->orientation;
}


std::string Voxel::toPovray(double ssize){
	std::stringstream out;

	out.precision(std::numeric_limits< double >::max_digits10);

	double *da = this->getPosition();
	double x1 = da[0], x2 = da[0] + ssize;
	double y1 = da[1], y2 = da[1] + ssize;

	out << "polygon {5, < " + std::to_string(x1) + ", " + std::to_string(y1) + ", 1.0>, "
					+ "< " + std::to_string(x1) + ", " + std::to_string(y2) + ", 1.0>, "
					+ "< " + std::to_string(x2) + ", " + std::to_string(y2) + ", 1.0>, "
					+ "< " + std::to_string(x2) + ", " + std::to_string(y1) + ", 1.0>, "
					+ "< " + std::to_string(x1) + ", " + std::to_string(y1) + ", 1.0>"
					+ " texture { pigment { color Black} } }";

	/*
			sRes += "\r\n  text { ttf \"timrom.ttf\" \"" + std::to_string(i) + "\" 1, 0 pigment { color White } scale 0.1 translate < ";
			for(unsigned char j=0; j<this->dimension; j++)
				sRes += std::to_string(da[j]+0.5*this->voxelSize) + ", ";
			sRes +=  "1.0003> }";
	*/

	return out.str();
}



std::string Voxel::toWolfram(double ssize, double asize){
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);

	double *da = this->getPosition();
	double x1 = da[0], x2 = da[0] + ssize;
	double y1 = da[1], y2 = da[1] + ssize;

	out << "Polygon[{ {" << x1 << ", " << y1 << "}, "
			<< "{" << x1 << ", " << y2 << "}, "
			<< "{" << x2 << ", " << y2 << "}, "
			<< "{" << x2 << ", " << y1 << "} }]";
	if (RSA_ANGULAR_DIMENSION > 0){
		out << "(* angles: [ " << this->getOrientation()[0] << ", " << (this->getOrientation()[0] + asize) << ") *)" << std::endl;
	}
	return out.str();
}



std::string Voxel::toString(){
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);

	out << " position: (";
	for (unsigned short i=0; i<RSA_SPATIAL_DIMENSION; i++){
		out << this->position[i];
		if (i<RSA_SPATIAL_DIMENSION-1)
			out << ", ";
	}
	out << ") orientation: (";
	for (unsigned short i=0; i<RSA_ANGULAR_DIMENSION; i++){
		out << this->orientation[i];
		if (i<RSA_ANGULAR_DIMENSION-1)
			out << ", ";
	}
	out << ")";
	return out.str();
}


void Voxel::store(std::ostream &f) const{
	unsigned short sd = RSA_SPATIAL_DIMENSION;
	unsigned short ad = RSA_ANGULAR_DIMENSION;
	f.write((char *)(&sd), sizeof(unsigned char));
	if (ad>0)
		f.write((char *)(&ad), sizeof(unsigned char));
	f.write((char *)(this->position), RSA_SPATIAL_DIMENSION*sizeof(double));
	if (ad>0)
		f.write((char *)(this->orientation.data()), RSA_ANGULAR_DIMENSION*sizeof(double));
}


void Voxel::restore(std::istream &f){
	unsigned char sd = RSA_SPATIAL_DIMENSION;
	unsigned char ad = RSA_ANGULAR_DIMENSION;

	f.read((char *)(&sd), sizeof(unsigned char));
	if (ad > 0)
		f.read((char *)(&ad), sizeof(unsigned char));

	if (sd!=RSA_SPATIAL_DIMENSION || ad!=RSA_ANGULAR_DIMENSION){
		std::cout << "[ERROR] cannot restore Voxel: incompatible dimensions: read " << f.gcount() << " bytes." << std::endl;
		return;
	}
	f.read((char *)(this->position), sd*sizeof(double));
	if (ad>0)
		f.read((char *)(this->orientation.data()), ad*sizeof(double));
}

