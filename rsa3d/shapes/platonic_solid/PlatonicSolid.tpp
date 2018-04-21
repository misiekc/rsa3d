//
// Created by PKua on 21.04.18.
//


template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::initClass(const std::string &attr) {
    SpecificSolid::calculateStatic(attr);
    Shape<3, 0>::setNeighbourListCellSize(2 * SpecificSolid::exsphereRadius);
    Shape<3, 0>::setVoxelSpatialSize(SpecificSolid::insphereRadius / M_SQRT2);
    Shape<3, 0>::setCreateShapeImpl([](RND *rnd) -> Shape<3, 0>* {
        return new SpecificSolid(Matrix<3, 3>::rotation(
                rnd->nextValue() * 2 * M_PI,
                std::asin(rnd->nextValue() * 2 - 1),
                rnd->nextValue() * 2 * M_PI));
    });
}

template<typename SpecificSolid>
int PlatonicSolid<SpecificSolid>::overlap(BoundaryConditions *bc, Shape<3, 0> *s) const {
    return 0;
}

template<typename SpecificSolid>
int PlatonicSolid<SpecificSolid>::pointInside(BoundaryConditions *bc, double *position,
                                              const std::array<double, 0> &orientation, double orientationRange) const {
    return 0;
}

template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::restore(std::istream &f) {
    Shape::restore(f);
}
template<typename SpecificSolid>
const Matrix<3, 3> &PlatonicSolid<SpecificSolid>::getOrientationMatrix() const {
    return orientation;
}
