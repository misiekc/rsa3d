//
// Created by PKua on 12.06.18.
//

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool ConvexShape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::pointInside(BoundaryConditions *bc, double *position) const{
    std::array<double, ANGULAR_DIMENSION> zeroAngle;
    zeroAngle.fill(0);
    return this->pointInside(bc, position, zeroAngle, 2*M_PI);
}

template<unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool ConvexShape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::voxelInside(BoundaryConditions *bc, const double *voxelPosition,
                                                                    const std::array<double, ANGULAR_DIMENSION> &orientation,
                                                                    double spatialSize, double angularSize) const {
    double position[SPATIAL_DIMENSION];
    int counterSize = 1 << SPATIAL_DIMENSION;
    bool isInside = true;
    for(int i=0; i<counterSize; i++){
        for(unsigned short j=0; j<SPATIAL_DIMENSION; j++){
            position[j] = voxelPosition[j] + Positioned<SPATIAL_DIMENSION>::offset[i][j]*spatialSize;
        }
        if(!this->pointInside(bc, position, orientation, angularSize)){
            isInside = false;
            break;
        }
    }
    return isInside;
}
