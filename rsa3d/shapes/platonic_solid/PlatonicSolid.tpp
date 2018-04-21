//
// Created by PKua on 21.04.18.
//

#include "../../Vector.h"

template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::initClass(const std::string &attr) {
    SpecificSolid::calculateStatic(attr);
    Shape::setNeighbourListCellSize(2 * SpecificSolid::exsphereRadius);
    Shape::setVoxelSpatialSize(SpecificSolid::insphereRadius / M_SQRT2);

    Shape::setCreateShapeImpl([](RND *rnd) -> Shape* {
        return new SpecificSolid(Matrix<3, 3>::rotation(
                rnd->nextValue() * 2 * M_PI,
                std::asin(rnd->nextValue() * 2 - 1),
                rnd->nextValue() * 2 * M_PI));
    });
}

template<typename SpecificSolid>
int PlatonicSolid<SpecificSolid>::overlap(BoundaryConditions *bc, Shape<3, 0> *s) const {
    /*auto *other = dynamic_cast<PlatonicSolid<SpecificSolid>*>(s);
    Vector<3> distance = vOtherPos - vThisPos;

    // Calculate edge Use std::array with appropriate size dictated by CRTP SpecificSolid
    using edge_array = decltype(SpecificSolid::edges);
    edge_array thisAxes, otherAxes;
    std::transform(SpecificSolid::edges.begin(), SpecificSolid::edges.end(), thisAxes,
                   [this](const Vector<3> &edge) {
                       return this->orientation * edge;
                   });
    std::transform(SpecificSolid::edges.begin(), SpecificSolid::edges.end(), otherAxes,
                   [other](const Vector<3> &edge) {
                       return other->orientation * edge;
                   });*/

    return 0;
}

template<typename SpecificSolid>
int PlatonicSolid<SpecificSolid>::pointInside(BoundaryConditions *bc, double *position,
                                              const std::array<double, 0> &orientation, double orientationRange) const {
    Vector<3> vThisPos(this->getPosition());
    Vector<3> vPointPos(position);

    return (vPointPos - vThisPos).norm2() <= SpecificSolid::insphereRadius;
}

template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::store(std::ostream &f) const {
    Shape::store(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            d = this->orientation(i, j);
            f.write((char *)(&d), sizeof(double));
        }
    }
}

template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::restore(std::istream &f) {
    Shape::restore(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            f.read((char *)&d, sizeof(double));
            this->orientation(i, j) = d;
        }
    }
}

template<typename SpecificSolid>
const Matrix<3, 3> &PlatonicSolid<SpecificSolid>::getOrientationMatrix() const {
    return orientation;
}

template<typename SpecificSolid>
bool PlatonicSolid<SpecificSolid>::isSeparatingAxis(const Vector<3> &axis, const PlatonicSolid &other,
                                                    const Vector<3> &distance) {
    double distanceProj = std::abs(distance * axis);
    return distanceProj > this->projectionHalfsize(axis) + other.projectionHalfsize(axis);
}