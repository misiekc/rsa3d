//
// Created by PKua on 21.04.18.
//

#include "../../Vector.h"

template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::initClass(const std::string &attr) {
    SpecificSolid::calculateStatic(attr);
    Shape::setNeighbourListCellSize(2 * SpecificSolid::circumsphereRadius);
    Shape::setVoxelSpatialSize(SpecificSolid::insphereRadius / M_SQRT2);

    Shape::setCreateShapeImpl([](RND *rnd) -> Shape* {
        return new SpecificSolid(Matrix<3, 3>::rotation(
                rnd->nextValue() * 2 * M_PI,
                std::asin(rnd->nextValue() * 2 - 1),
                rnd->nextValue() * 2 * M_PI));
    });
}

template<typename SpecificSolid>
template<size_t SIZE>
std::array<Vector<3>, SIZE> PlatonicSolid<SpecificSolid>::applyOrientation(const std::array<Vector<3>, SIZE> &vectors) const {
    std::array<Vector<3>, SIZE> result;
    std::transform(vectors.begin(), vectors.end(), result.begin(), [this](const Vector<3> &v) {
        return this->orientation * v;
    });
    return result;
}

template<typename SpecificSolid>
int PlatonicSolid<SpecificSolid>::overlap(BoundaryConditions *bc, Shape<3, 0> *s) const {
    SpecificSolid other = dynamic_cast<SpecificSolid&>(*s);   // Make a copy
    this->applyBC(bc, &other);

    // TODO maybe store rotated axes in SpecificSolid instances?
    auto thisFaceAxes = this->applyOrientation(SpecificSolid::orientedFaceAxes);
    auto thisEdgeAxes = this->applyOrientation(SpecificSolid::orientedEdgeAxes);
    auto otherFaceAxes = other.applyOrientation(SpecificSolid::orientedFaceAxes);
    auto otherEdgeAxes = other.applyOrientation(SpecificSolid::orientedEdgeAxes);

    Vector<3> distance = Vector<3>(other.getPosition()) - Vector<3>(this->getPosition());
    auto thisSpecific = static_cast<const SpecificSolid *>(this);
    // Face axes for this
    for (const auto &face : thisFaceAxes)
        if (thisSpecific->isSeparatingAxis(face, other, distance))
            return 0;
    // Face axes for other
    for (const auto &face : otherFaceAxes)
        if (thisSpecific->isSeparatingAxis(face, other, distance))
            return 0;
    // Egde axes cross products
    for (const auto &thisEdge : thisEdgeAxes)
        for (const auto &otherEdge : otherEdgeAxes)
            if (thisSpecific->isSeparatingAxis(thisEdge ^ otherEdge, other, distance))
                return 0;

    return 1;
}

template<typename SpecificSolid>
int PlatonicSolid<SpecificSolid>::pointInside(BoundaryConditions *bc, double *position,
                                              const std::array<double, 0> &orientation, double orientationRange) const {
    Vector<3> vThisPos(this->getPosition());
    Vector<3> vPointPos(position);

    return (vPointPos - vThisPos).norm2() <= 4 * std::pow(SpecificSolid::insphereRadius, 2);
}

template<typename SpecificSolid>
bool PlatonicSolid<SpecificSolid>::isSeparatingAxis(const Vector<3> &axis, const SpecificSolid &other,
                                                    const Vector<3> &distance) const {
    auto &thisSpecific = static_cast<const SpecificSolid&>(*this);

    double distanceProj = std::abs(distance * axis);
    return distanceProj > thisSpecific.projectionHalfsize(axis) + other.projectionHalfsize(axis);
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
