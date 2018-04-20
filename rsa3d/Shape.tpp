/*
 * Shape.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include <iostream>

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::voxelSpatialSize = 0;

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::voxelAngularSize = 2*M_PI;

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::neighbourListCellSize = 0;

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
typename Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::create_shape_fun_ptr
        Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::createShapeImpl = nullptr;


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setVoxelSpatialSize(double size) {
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::voxelSpatialSize = size;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setVoxelAngularSize(double size) {
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::voxelAngularSize = size;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setNeighbourListCellSize(double size) {
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::neighbourListCellSize = size;
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVoxelSpatialSize() {
	return Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::voxelSpatialSize;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVoxelAngularSize() {
	return Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::voxelAngularSize;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getNeighbourListCellSize() {
	return Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::neighbourListCellSize;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
const typename Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::create_shape_fun_ptr
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getCreateShapeImpl() {
    return createShapeImpl;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::Shape() : Positioned<SPATIAL_DIMENSION>(){
    this->orientation.fill(0);
	this->no = 0;
	this->time = 0.0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVolume() const {
    return 1.;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::array<double, ANGULAR_DIMENSION> Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getOrientation() const{
    return this->orientation;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setOrientation(const double *orientation) {
    std::copy(orientation, orientation + ANGULAR_DIMENSION, this->orientation.begin());
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::rotate(double* v){
    std::array<double, ANGULAR_DIMENSION> orientation;
	for(unsigned short i=0; i<ANGULAR_DIMENSION; i++)
		orientation[i] = this->getOrientation()[i] + v[i];
    this->setOrientation(orientation.data());
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::pointInside(BoundaryConditions *bc, double *position) const{
    std::array<double, ANGULAR_DIMENSION> zeroAngle;
    zeroAngle.fill(0);
	return this->pointInside(bc, position, zeroAngle, 2*M_PI);
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::minDistance(Shape *s) const{
	return 0.0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toString() const{
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
    // TODO move to Positioned::store?
	f.write((char *)(this->getPosition()), SPATIAL_DIMENSION*sizeof(double));
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
        throw std::runtime_error(
                std::string("[ERROR] cannot restore Shape: incompatible dimensions: read ")
                + std::to_string(f.gcount())
                + " bytes.");
	}
	f.read((char *)(&this->no), sizeof(int));
	f.read((char *)(&this->time), sizeof(double));

	// TODO move to Positioned::restore?
	double position[SPATIAL_DIMENSION];
	f.read((char *)(position), sd*sizeof(double));
	this->setPosition(position);

	if (ad>0)
		f.read((char *)(this->orientation.data()), ad*sizeof(double));
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::applyBC(BoundaryConditions *bc,
                                                          Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *second) const {
    double translation[SPATIAL_DIMENSION];
    bc->getTranslation(translation, this->getPosition(), second->getPosition());
    second->translate(translation);
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setCreateShapeImpl(Shape *(*const fptr)(RND *)) {
    createShapeImpl = fptr;
}