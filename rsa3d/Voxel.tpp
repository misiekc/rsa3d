/*
 * Voxel.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */
#include <sstream>
#include <stdexcept>

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::Voxel(){
	this->index = 0;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
	this->depth = 0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::Voxel(double* pos, const std::array<double, ANGULAR_DIMENSION> &angle){
	std::copy(pos, pos+SPATIAL_DIMENSION, this->position);
	this->orientation = angle;
	this->index = 0;
	this->missCounter = 0;
	this->lastAnalyzed = 0;
	this->depth = 0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::miss(){
	this->missCounter++;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getMissCounter(){
	return this->missCounter;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::resetMissCounter(){
	this->missCounter = 0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::isInside(double *pos, double size){
	for(int i=0; i<SPATIAL_DIMENSION; i++){
		if (pos[i]<this->position[i])
			return false;
		if (pos[i]>=(this->position[i]+size))
			return false;
	}
	return true;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::isInside(double *pos, double size, const std::array<double, ANGULAR_DIMENSION> &angle, double asize){
	for(int i=0; i<SPATIAL_DIMENSION; i++){
		if (pos[i]<this->position[i])
			return false;
		if (pos[i]>=(this->position[i]+size))
			return false;
	}

	for(int i=0; i<ANGULAR_DIMENSION; i++){
		if (angle[i]<this->orientation[i])
			return false;
		if (angle[i]>=(this->orientation[i]+asize))
			return false;
	}

	return true;
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double* Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getPosition(){
	return this->position;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::array<double, ANGULAR_DIMENSION> Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getOrientation(){
	return this->orientation;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toPovray(double ssize){
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


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toWolfram(double ssize, double asize){
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);

	double *da = this->getPosition();
	double x1 = da[0], x2 = da[0] + ssize;
	double y1 = da[1], y2 = da[1] + ssize;

	out << "Polygon[{ {" << x1 << ", " << y1 << "}, "
			<< "{" << x1 << ", " << y2 << "}, "
			<< "{" << x2 << ", " << y2 << "}, "
			<< "{" << x2 << ", " << y1 << "} }]";
	if (ANGULAR_DIMENSION > 0){
		out << "(* angles: [ " << this->getOrientation()[0] << ", " << (this->getOrientation()[0] + asize) << ") *)" << std::endl;
	}
	return out.str();
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toString(){
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);

	out << "index: " << this->index << " position: (";
	for (unsigned short i=0; i<SPATIAL_DIMENSION; i++){
		out << this->position[i];
		if (i<SPATIAL_DIMENSION-1)
			out << ", ";
	}
	out << ") orientation: (";
	for (unsigned short i=0; i<ANGULAR_DIMENSION; i++){
		out << this->orientation[i];
		if (i<ANGULAR_DIMENSION-1)
			out << ", ";
	}
	out << ")";
	return out.str();
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::store(std::ostream &f) const{
	unsigned short sd = SPATIAL_DIMENSION;
	unsigned short ad = ANGULAR_DIMENSION;
	f.write((char *)(&sd), sizeof(unsigned char));
	if (ad>0)
		f.write((char *)(&ad), sizeof(unsigned char));
	f.write((char *)(this->position), SPATIAL_DIMENSION*sizeof(double));
	if (ad>0)
		f.write((char *)(this->orientation), ANGULAR_DIMENSION*sizeof(double));
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::restore(std::istream &f){
	unsigned char sd = SPATIAL_DIMENSION;
	unsigned char ad = ANGULAR_DIMENSION;

	f.read((char *)(&sd), sizeof(unsigned char));
	if (ad > 0)
		f.read((char *)(&ad), sizeof(unsigned char));

	if (sd!=SPATIAL_DIMENSION || ad!=ANGULAR_DIMENSION){
		std::cout << "[ERROR] cannot restore Voxel: incompatible dimensions: read " << f.gcount() << " bytes." << std::endl;
		return;
	}
	f.read((char *)(this->position), sd*sizeof(double));
	if (ad>0)
		f.read((char *)(this->orientation.data()), ad*sizeof(double));
}

