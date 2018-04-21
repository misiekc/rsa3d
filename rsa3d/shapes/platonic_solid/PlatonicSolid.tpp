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

    // Calculate edge and face axes for current positions. Use std::array with appropriate size dictated by CRTP
    // SpecificSolid
    // TODO maybe store rotated axes in SpecificSolid instances?
    using face_array = decltype(SpecificSolid::faceAxes);
    using edge_array = decltype(SpecificSolid::edgeAxes);
    face_array thisFaceAxes, otherFaceAxes;
    edge_array thisEdgeAxes, otherEdgeAxes;

    auto this_rotater = [this](const Vector<3> &edge) {
        return this->orientation * edge;
    };

    auto other = dynamic_cast<PlatonicSolid<SpecificSolid>*>(s);
    auto other_rotater = [other](const Vector<3> &edge) {
        return other->orientation * edge;
    };

    std::transform(SpecificSolid::faceAxes.begin(), SpecificSolid::faceAxes.end(), thisFaceAxes.begin(), this_rotater);
    std::transform(SpecificSolid::faceAxes.begin(), SpecificSolid::faceAxes.end(), thisEdgeAxes.begin(), this_rotater);
    std::transform(SpecificSolid::edgeAxes.begin(), SpecificSolid::edgeAxes.end(), otherFaceAxes.begin(), other_rotater);
    std::transform(SpecificSolid::edgeAxes.begin(), SpecificSolid::edgeAxes.end(), otherEdgeAxes.begin(), other_rotater);

    auto thisSpecific = static_cast<const SpecificSolid *>(this);
    Vector<3> distance = Vector<3>(s->getPosition()) - Vector<3>(this->getPosition());
    // Face axes for this
    for (const auto &face : thisFaceAxes)
        if (thisSpecific->isSeparatingAxis(face, *other, distance))
            return 0;
    // Face axes for other
    for (const auto &face : otherFaceAxes)
        if (thisSpecific->isSeparatingAxis(face, *other, distance))
            return 0;
    // Egde axes cross products
    for (const auto &thisEdge : thisEdgeAxes)
        for (const auto &otherEdge : otherEdgeAxes)
            if (thisSpecific->isSeparatingAxis(thisEdge ^ otherEdge, *other, distance))
                return 0;

    return 1;
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
                                                    const Vector<3> &distance) const {
    auto &thisSpecific = static_cast<const SpecificSolid&>(*this);
    auto &otherSpecific = static_cast<const SpecificSolid&>(other);

    double distanceProj = std::abs(distance * axis);
    return distanceProj > thisSpecific.projectionHalfsize(axis) + otherSpecific.projectionHalfsize(axis);
}