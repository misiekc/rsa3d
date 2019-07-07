/*
 * Shape.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include <iostream>
#include "Shape.h"


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION> Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::shapeStaticInfo;

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::earlyRejectionEnabled = true;


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::Shape() : Positioned<SPATIAL_DIMENSION>(){
    this->orientation.fill(0);
	this->no = 0;
	this->time = 0.0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
const Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::overlap(BoundaryConditions<SPATIAL_DIMENSION> *bc,
						std::vector<const Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *shapes) const{

	for(const Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s: *shapes)
		if (this->overlap(bc, s))
			return s;

	return nullptr;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getVolume(unsigned short dim) const {
    // By default, particles have unity volume and do not support packings in dimensions other than their own
    if (dim != SPATIAL_DIMENSION)
        throw std::runtime_error (" supports only " + std::to_string(SPATIAL_DIMENSION) + "D packings");
    return 1.;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
Orientation<ANGULAR_DIMENSION> Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getOrientation() const{
    return this->orientation;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setOrientation(
        const Orientation<ANGULAR_DIMENSION> &orientation) {
	this->orientation = orientation;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::rotate(const Orientation<ANGULAR_DIMENSION> &v){
    Orientation<ANGULAR_DIMENSION> orientation;
	for(unsigned short i=0; i<ANGULAR_DIMENSION; i++)
		orientation[i] = this->getOrientation()[i] + v[i];
    this->setOrientation(orientation);
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::minDistance(const Shape *s) const{
	return 0.0;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::string Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toString() const{
	return std::to_string(this->no);
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
    double position[SPATIAL_DIMENSION];
    this->getPosition().copyToArray(position);
	f.write((char *)(position), SPATIAL_DIMENSION*sizeof(double));
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
	this->setPosition(Vector<SPATIAL_DIMENSION>(position));

	if (ad>0)
		f.read((char *)(this->orientation.data()), ad*sizeof(double));
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::applyBC(BoundaryConditions<SPATIAL_DIMENSION> *bc,
                                                          Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *second) const {
    second->translate(bc->getTranslation(this->getPosition(), second->getPosition()));
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setShapeStaticInfo(
        ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION> shapeStaticInfo) {
    shapeStaticInfo.throwIfIncomplete();
    Shape::shapeStaticInfo = std::move(shapeStaticInfo);
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
typename Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::EarlyRejectionResult
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::overlapEarlyRejection(
        BoundaryConditions<SPATIAL_DIMENSION> *bc, const Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s) const {
    if (!earlyRejectionEnabled)
        return UNKNOWN;

    double distance2 = bc->distance2(this->getPosition(), s->getPosition());

    if (distance2 < 4*std::pow(shapeStaticInfo.getInsphereRadius(), 2))
        return TRUE;
    else if (distance2 > 4*std::pow(shapeStaticInfo.getCircumsphereRadius(), 2))
        return FALSE;
    else
        return UNKNOWN;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setEarlyRejectionEnabled(bool earlyRejectionEnabled) {
    Shape::earlyRejectionEnabled = earlyRejectionEnabled;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
typename Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::EarlyRejectionResult
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::voxelInsideEarlyRejection(
        BoundaryConditions<SPATIAL_DIMENSION> *bc, const Vector<SPATIAL_DIMENSION> &voxelPosition,
        const Orientation<ANGULAR_DIMENSION> &orientation, double spatialSize, double angularSize) const {
    if (!earlyRejectionEnabled)
        return UNKNOWN;

    double halfSpatialSize = 0.5*spatialSize;
    auto voxelTranslation = bc->getTranslation(this->getPosition(), voxelPosition);
    auto spatialCenter = voxelPosition + voxelTranslation + Vector<SPATIAL_DIMENSION>(halfSpatialSize);
    auto relativeSpatialCenter = spatialCenter - this->getPosition();

    /* Both early rejection mechanisms */
    /* EarlyRejectionResult result = sphereOnVoxelEarlyRejection(relativeSpatialCenter, halfSpatialSize);
    if (result != UNKNOWN)
        return result;
    else
        return furthestCornerEarlyRejection(relativeSpatialCenter, halfSpatialSize);*/

    /* Sphere on voxel early rejection only */
    /* return sphereOnVoxelEarlyRejection(relativeSpatialCenter, halfSpatialSize); */

    /* Furthest corner early rejection only - seems the best option */
    return furthestCornerEarlyRejection(relativeSpatialCenter, halfSpatialSize);
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
typename Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::EarlyRejectionResult
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::sphereOnVoxelEarlyRejection(
        const Vector<SPATIAL_DIMENSION> &relativeSpatialCenter, double halfSpatialSize) const {
    double spatialCenterDistance2 = relativeSpatialCenter.norm2();
    double voxelCircumsphereRadius = halfSpatialSize * std::sqrt(SPATIAL_DIMENSION);
    double voxelInsphereRadius = halfSpatialSize;

    if (spatialCenterDistance2 < std::pow(shapeStaticInfo.getExclusionZoneMinSpan() - voxelCircumsphereRadius, 2))
        return TRUE;
    else if (spatialCenterDistance2 > std::pow(shapeStaticInfo.getExclusionZoneMaxSpan() - voxelInsphereRadius, 2))
        return FALSE;
    else
        return UNKNOWN;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
typename Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::EarlyRejectionResult
Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::furthestCornerEarlyRejection(
        const Vector<SPATIAL_DIMENSION> &relativeSpatialCenter, double halfSpatialSize) const {
    auto furthestVoxelCorner = relativeSpatialCenter;
    for (std::size_t j = 0; j < SPATIAL_DIMENSION; j++) {
        if (relativeSpatialCenter[j] > 0)
            furthestVoxelCorner[j] += halfSpatialSize;
        else
            furthestVoxelCorner[j] -= halfSpatialSize;
    }

    double furthestVoxelCornerDistance2 = furthestVoxelCorner.norm2();
    if (furthestVoxelCornerDistance2 < std::pow(shapeStaticInfo.getExclusionZoneMinSpan(), 2))
        return TRUE;
    else if (furthestVoxelCornerDistance2 > std::pow(shapeStaticInfo.getExclusionZoneMaxSpan(), 2))
        return FALSE;
    else
        return UNKNOWN;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setCircumsphereRadius(double circumsphereRadius) {
    Expects(circumsphereRadius > 0);
    Expects(circumsphereRadius < std::numeric_limits<double>::infinity());
    this->circumsphereRadius = circumsphereRadius;
    if (this->neighbourListCellSize == NOT_SPECIFIED)
        this->neighbourListCellSize = 2*circumsphereRadius;
    if (this->exclusionZoneMaxSpan == NOT_SPECIFIED)
        this->exclusionZoneMaxSpan = 2*circumsphereRadius;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setInsphereRadius(double insphereRadius) {
    Expects(insphereRadius > 0);
    Expects(insphereRadius < std::numeric_limits<double>::infinity());
    this->insphereRadius = insphereRadius;
    if (this->voxelSpatialSize == NOT_SPECIFIED)
        this->voxelSpatialSize = 2*insphereRadius/std::sqrt(SPATIAL_DIMENSION);
    if (this->exclusionZoneMinSpan == NOT_SPECIFIED)
        this->exclusionZoneMinSpan = 2*insphereRadius;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setNeighbourListCellSize(double neighbourListCellSize) {
    Expects(neighbourListCellSize > 0);
    Expects(neighbourListCellSize < std::numeric_limits<double>::infinity());
    ShapeStaticInfo::neighbourListCellSize = neighbourListCellSize;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setVoxelSpatialSize(double voxelSpatialSize) {
    Expects(voxelSpatialSize > 0);
    Expects(voxelSpatialSize < std::numeric_limits<double>::infinity());
    ShapeStaticInfo::voxelSpatialSize = voxelSpatialSize;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setVoxelAngularSize(double voxelAngularSize) {
    Expects(voxelAngularSize > 0);
    Expects(voxelAngularSize <= 2*M_PI);
    ShapeStaticInfo::voxelAngularSize = voxelAngularSize;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setExclusionZoneMinSpan(double exclusionZoneMinSpan) {
    Expects(exclusionZoneMinSpan > 0);
    Expects(exclusionZoneMinSpan < std::numeric_limits<double>::infinity());
    ShapeStaticInfo::exclusionZoneMinSpan = exclusionZoneMinSpan;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setExclusionZoneMaxSpan(double exclusionZoneMaxSpan) {
    Expects(exclusionZoneMaxSpan > 0);
    Expects(exclusionZoneMaxSpan < std::numeric_limits<double>::infinity());
    ShapeStaticInfo::exclusionZoneMaxSpan = exclusionZoneMaxSpan;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setSupportsSaturation(bool supportsSaturation) {
    ShapeStaticInfo::supportsSaturation = supportsSaturation;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setCreateShapeImpl(
        ShapeStaticInfo::create_shape_fun_ptr createShapeImpl) {
    ShapeStaticInfo::createShapeImpl = createShapeImpl;
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
template<typename ConcreteShape>
void ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::setDefaultCreateShapeImpl() {
    this->createShapeImpl = [](RND *rnd) {
        return static_cast<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *>(new ConcreteShape());
    };
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void ShapeStaticInfo<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::throwIfIncomplete() const {
    if (this->circumsphereRadius == NOT_SPECIFIED)
        throw std::runtime_error("Circumsphere radius has not been set!");
    else if (this->insphereRadius == NOT_SPECIFIED)
        throw std::runtime_error("Insphere radius has not been set!");
    else if (this->createShapeImpl == nullptr)
        throw std::runtime_error("Create shape function implementation radius has not been set!");
}
