//
// Created by ciesla on 12/23/22.
//

#include "NewDiscreteAngleVoxelList.h"
#include "utils/OMPMacros.h"

#include <cmath>


NewDiscreteAngleVoxelList::NewDiscreteAngleVoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange, double requestedAngularVoxelSize, const std::vector<Orientation<1>> &orientations) : VoxelList(dim, packingSpatialSize, requestedSpatialVoxelSize, shapeAngularRange, requestedAngularVoxelSize){
    this->allowedOrientations = orientations;
    this->voxelMap = new double[1];
    this->u01Distribution = new std::uniform_real_distribution<double>(0.0, 1.0);
    this->orientationsVectorRange = this->angularRange;
    this->orientationsVector = new std::vector<RSAOrientation >[1];
    this->orientationsVector[0] = this->allowedOrientations;
    this->orientationsHalfVector = nullptr;
 }

NewDiscreteAngleVoxelList::~NewDiscreteAngleVoxelList() {
    delete this->voxelMap;
    delete this->u01Distribution;
    if (this->orientationsVector!= nullptr)
        delete[] this->orientationsVector;
    if (this->orientationsHalfVector!= nullptr)
        delete[] this->orientationsHalfVector;

}

void NewDiscreteAngleVoxelList::allocateVoxels(size_t size){
    delete[] this->voxels;

    this->voxels = new Voxel*[size];
    for(size_t i=0; i<size; i++){
        this->voxels[i] = nullptr;
    }
}

std::vector<RSAOrientation> NewDiscreteAngleVoxelList::getPossibleOrientations(double start, double range) const {
    std::vector<RSAOrientation> v;

    size_t i = start / range;
    if (range == this->orientationsVectorRange){
        v = this->orientationsVector[i];
    }else if (range == this->orientationsHalfVectorRange){
        v = this->orientationsHalfVector[i];
    }else {
        double end = start + range;
        for (const RSAOrientation orientation: this->allowedOrientations) {
            if (orientation[0] >= start && orientation[0] < end)
                v.push_back(orientation);
        }
    }
    return v;
}

double NewDiscreteAngleVoxelList::getAngularVoxelSize(size_t i) const{
    if (!this->voxelsInitialized)
        return this->angularRange;
    else{
        std::vector<RSAOrientation> v = this->getPossibleOrientations(this->voxels[i]->getOrientation()[0], this->angularVoxelSize);
        return this->angularRange * static_cast<double>(v.size()) / static_cast<double>(this->allowedOrientations.size());
    }
}

double NewDiscreteAngleVoxelList::getVoxelVolume(size_t i) const{
    double spatialVolume = std::pow(this->getSpatialVoxelSize(), this->surfaceDimension);
    double angularVolume = std::pow(this->getAngularVoxelSize(i), RSA_ANGULAR_DIMENSION);
    return (spatialVolume*angularVolume);
}

void NewDiscreteAngleVoxelList::createVoxelMap() {
    double sum = 0.0;
    delete this->voxelMap;
    this->voxelMap = new double[this->length];
    for (size_t i = 0; i < this->length; i++) {
        sum += this->getVoxelVolume(i);
        this->voxelMap[i] = sum;
    }
    for (size_t i = 0; i < this->length; i++) {
        this->voxelMap[i] /= sum;
    }
}

void NewDiscreteAngleVoxelList::createOrientationsMap(double range) {

    if (this->orientationsVector != nullptr)
        delete[] this->orientationsVector;

    if (this->orientationsHalfVector != nullptr && this->orientationsHalfVectorRange == range) {
        this->orientationsVector = this->orientationsHalfVector;
        this->orientationsHalfVector = nullptr;
    } else {
        size_t iMax = this->angularRange / range + 1;
        this->orientationsVector = new std::vector<RSAOrientation>[iMax];
        for (const RSAOrientation orientation: this->allowedOrientations) {
            size_t i = orientation[0] / range;
            this->orientationsVector[i].push_back(orientation);
        }
    }
    this->orientationsVectorRange = range;

    if (this->orientationsHalfVector != nullptr)
        delete[] this->orientationsHalfVector;
    size_t iMax = this->angularRange / ((range / 2.0) * this->dxFactor) + 1;
    this->orientationsHalfVector = new std::vector<RSAOrientation>[iMax];
    for (const RSAOrientation orientation: this->allowedOrientations) {
        size_t i = orientation[0] / ((range / 2.0) * this->dxFactor);
        this->orientationsHalfVector[i].push_back(orientation);
    }
    this->orientationsHalfVectorRange = ((range / 2.0) * this->dxFactor);
}

void NewDiscreteAngleVoxelList::getRandomEntry(RSAVector *position, RSAOrientation *orientation, Voxel **v, RND *rnd) {
    double d = rnd->nextValue();
    size_t i1=0, i2=this->length-1, i;
    while(i2>i1+1 && i2>0){
        i = 0.5*(i1 + i2);
        if (this->voxelMap[i] > d)
            i2 = i-1;
        else
            i1 = i;
    }
    *v = this->voxels[i1];
    this->getRandomPositionAndOrientation(position, orientation, *v, rnd);
}

void NewDiscreteAngleVoxelList::getRandomPositionAndOrientation(RSAVector *position, RSAOrientation *orientation, Voxel *v, RND *rnd){
    RSAVector vpos = v->getPosition();

    for (unsigned short i=0; i < this->surfaceDimension; i++)
        (*position)[i] = vpos[i] + rnd->nextValue(this->u01Distribution)*this->spatialVoxelSize;
    for (unsigned short i=this->surfaceDimension; i < RSA_SPATIAL_DIMENSION; i++)
        (*position)[i] = 0.0;

    std::vector<RSAOrientation> vo = this->getPossibleOrientations(v->getOrientation()[0], this->angularVoxelSize);
    RSAOrientation o = vo[static_cast<size_t>(vo.size()*rnd->nextValue(this->u01Distribution))];
    for (unsigned short i=0; i < RSA_ANGULAR_DIMENSION; i++)
        (*orientation)[i] = o[i];
}


unsigned short NewDiscreteAngleVoxelList::splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc){
    if (!this->voxelsInitialized){
        this->createOrientationsMap(this->initialAngularVoxelSize);
        this->initVoxels(bc, nl);
        this->createVoxelMap();
        return VoxelList::NO_SPLIT_BUT_INITIALIZED;
    }
    unsigned short uRet = VoxelList::splitVoxels(minDx, maxVoxels, nl, bc);
    this->createVoxelMap();
    this->createOrientationsMap(this->angularVoxelSize);
    return uRet;
}

double NewDiscreteAngleVoxelList::getVoxelsVolume() const{
    double result = 0;
    for(size_t i = 0; i< this->length; i++){
        double s = 1.0;
        Voxel *v = this->voxels[i];
        RSAVector position = v->getPosition();
        for(unsigned short j=0; j<this->surfaceDimension; j++){
            if (position[j]+this->spatialVoxelSize > this->spatialRange){
                s *= this->spatialRange - position[j];
            }else{
                s *= this->spatialVoxelSize;
            }
        }
        double a = 1.0;
        for(unsigned char j=0; j<RSA_ANGULAR_DIMENSION; j++){
            a *= this->getAngularVoxelSize(i);
            a /= this->angularRange;
        }
        result += s*a;
    }
    return result;
}

bool NewDiscreteAngleVoxelList::analyzeVoxel(Voxel *v, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc, double spatialSize, double angularSize, unsigned short depth) const{
    if (v->getOrientation()[0]>this->angularRange)
        return true;
    std::vector<RSAOrientation> vo = this->getPossibleOrientations(v->getOrientation()[0], angularSize);
    if (vo.size()==0)
        return true;
    else
        return VoxelList::analyzeVoxel(v, nl, bc, spatialSize, angularSize, depth);
}


