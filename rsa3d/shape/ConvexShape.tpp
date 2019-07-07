//
// Created by PKua on 12.06.18.
//

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool ConvexShape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::pointInside(BoundaryConditions<SPATIAL_DIMENSION> *bc,
                                                                    const Vector<SPATIAL_DIMENSION> &position) const{
    Orientation<ANGULAR_DIMENSION> zeroAngle;
    zeroAngle.fill(0);
    return this->pointInside(bc, position, zeroAngle, 2*M_PI);
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool ConvexShape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::voxelInside(BoundaryConditions<SPATIAL_DIMENSION> *bc,
                                                                    const Vector<SPATIAL_DIMENSION> &voxelPosition,
                                                                    const Orientation<ANGULAR_DIMENSION> &orientation,
                                                                    double spatialSize, double angularSize) const {
    switch(this->voxelInsideEarlyRejection(bc, voxelPosition, orientation, spatialSize, angularSize)) {
        case EarlyRejectionResult::TRUE:      return true;
        case EarlyRejectionResult::FALSE:     return false;
        case EarlyRejectionResult::UNKNOWN:   break;
    }

    int counterSize = 1 << SPATIAL_DIMENSION;
    bool isInside = true;
    for(int i=0; i<counterSize; i++){
        Vector<SPATIAL_DIMENSION> position;
        for(unsigned short j=0; j<SPATIAL_DIMENSION; j++){
            position[j] = voxelPosition[j] + Positioned<SPATIAL_DIMENSION>::offset[i][j]*spatialSize;
        }
        if(!this->pointInside(bc, position, orientation, angularSize)){
            isInside = false;
            break;
        }
    }

    /* Temporary consistency check. Will be deleting after proper test implementation */
    /*auto result = this->voxelInsideEarlyRejection(bc, voxelPosition, orientation, spatialSize, angularSize);
    if (result != EarlyRejectionResult::UNKNOWN) {
        bool boolResult = (result == EarlyRejectionResult::TRUE);

        if (boolResult != isInside) {
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << "Early rejection mismatch!" << std::endl;
            std::cout << "Early: " << boolResult << ", real: " << isInside << std::endl;
            std::cout << "V position: " << voxelPosition << ", spatialSize: " << spatialSize << std::endl;
            std::cout << "Circumpshere: " << Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getCircumsphereRadius();
            std::cout << ", insphere: " << Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getInsphereRadius();
            std::cout << "shape:" << std::endl << this->toWolfram() << std::endl;
        }
    }*/

    return isInside;
}
