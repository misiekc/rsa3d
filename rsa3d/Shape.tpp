/*
 * Shape.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include <iostream>

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVoxelAngularSize(){
	return 2*M_PI;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::Shape() : Positioned<SPATIAL_DIMENSION>(){
    this->orientation.fill(0);
	this->no = 0;
	this->time = 0.0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::Shape(const Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> & other) :
        Positioned<SPATIAL_DIMENSION>(other), orientation(other.orientation) {
	this->no = 0;
	this->time = 0.0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::~Shape() {

}

// Version for shapes with non-zero angular dimension
/*template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
template <unsigned short AD>
typename std::enable_if<AD != 0, const double*>::type
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getOrientation() const{
    return this->orientation;
}

// Version for shapes with no angular dimension
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
template <unsigned short AD>
typename std::enable_if<AD == 0, const double*>::type
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getOrientation() const{
    return nullptr;     // return nullptr so dereferencing will fail
}*/

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
const double* Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getOrientation() const{
    return this->orientation.data();
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setOrientation(const double *orientation) {
    std::copy(orientation, orientation + ANGULAR_DIMENSION, this->orientation.begin());
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::translate(double* v){
	for(unsigned short i=0; i<SPATIAL_DIMENSION; i++){
		this->position[i] += v[i];
	}
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::rotate(double* v){
	for(unsigned short i=0; i<ANGULAR_DIMENSION; i++){
		this->orientation[i] += v[i];
	}
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::pointInside(BoundaryConditions *bc, double *position){
	return this->pointInside(bc, position, 0, 2*M_PI);
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::minDistance(Shape *s){
	return 0.0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toString(){
	return "";
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toPovray() const{
	return "";
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toWolfram() const{
	return "";
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::store(std::ostream &f) const{
	unsigned short sd = SPATIAL_DIMENSION;
	unsigned short ad = ANGULAR_DIMENSION;
	f.write((char *)(&sd), sizeof(unsigned char));
	if (ad>0)
		f.write((char *)(&ad), sizeof(unsigned char));
	f.write((char *)(&this->no), sizeof(int));
	f.write((char *)(&this->time), sizeof(double));
	f.write((char *)(this->position), SPATIAL_DIMENSION*sizeof(double));
	if (ad>0)
		f.write((char *)(this->orientation.data()), ANGULAR_DIMENSION*sizeof(double));
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::restore(std::istream &f){
	unsigned char sd = SPATIAL_DIMENSION;
	unsigned char ad = ANGULAR_DIMENSION;

	f.read((char *)(&sd), sizeof(unsigned char));
	if (f.gcount()==0){ // end of file
		return;
	}
	if (ad > 0)
		f.read((char *)(&ad), sizeof(unsigned char));

	if (sd!=SPATIAL_DIMENSION || ad!=ANGULAR_DIMENSION){
		std::cout << "[ERROR] cannot restore Shape: incompatible dimensions: read " << f.gcount() << " bytes." << std::endl;
		return;
	}
	f.read((char *)(&this->no), sizeof(int));
	f.read((char *)(&this->time), sizeof(double));
	f.read((char *)(this->position), sd*sizeof(double));
	if (ad>0)
		f.read((char *)(this->orientation.data()), ad*sizeof(double));
}

/*
template <unsigned short SPATIAL_DIMENSION>
void Shape<SPATIAL_DIMENSION>::vectorTranslate(const Vector<SPATIAL_DIMENSION> &translation) {
	double arr[SPATIAL_DIMENSION];
	translation.copyToArray(arr);
	translate(arr);
}
*/

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::applyBC(BoundaryConditions *bc, Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *second) {
    double translation[SPATIAL_DIMENSION];
    bc->getTranslation(translation, this->getPosition(), second->getPosition());
    second->translate(translation);
}

/*
template <unsigned short SPATIAL_DIMENSION>
double* Shape<SPATIAL_DIMENSION>::applyBC(double *res, BoundaryConditions *bc, double *pointToTranslate) {
	double* Shape<SPATIAL_DIMENSION>::applyBC(double *res, BoundaryConditions *bc, double *pointToTranslate) {
	bc->getTranslation(res, this->position, pointToTranslate);
	return res;
}
*/
