//
// Created by PKua on 21.04.18.
//

#include "../../Vector.h"
#include "SATOverlap.h"
#include "TriTriOverlap.h"
#include <sstream>

template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::initClass(const std::string &attr) {
    SpecificSolid::calculateStatic(attr);
    normalizeAxes();

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
void PlatonicSolid<SpecificSolid>::normalizeAxes() {
    auto normalizer = [](const Vector<3> &v) {
        return v / v.norm();
    };
    std::transform(SpecificSolid::orientedEdgeAxes.begin(), SpecificSolid::orientedEdgeAxes.end(),
                   SpecificSolid::orientedEdgeAxes.begin(), normalizer);
    std::transform(SpecificSolid::orientedFaceAxes.begin(), SpecificSolid::orientedFaceAxes.end(),
                   SpecificSolid::orientedFaceAxes.begin(), normalizer);
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
    SpecificSolid other = dynamic_cast<SpecificSolid&>(*s);     // Make a copy
    this->applyBC(bc, &other);
    SATOverlap<SpecificSolid> satOverlap;
    return satOverlap.overlap(this, &other);
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
    Matrix<3, 3> orientation;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            f.read((char *)&d, sizeof(double));
            orientation(i, j) = d;
        }
    }
    this->setOrientationMatrix(orientation);
}

template<typename SpecificSolid>
const Matrix<3, 3> &PlatonicSolid<SpecificSolid>::getOrientationMatrix() const { return orientation; }

template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::setOrientationMatrix(const Matrix<3, 3> &orientation) {
    this->orientation = orientation;
    calculateVerticesAndAxes();
}

template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::setPosition(const double *position) {
    Vector<3> oldPos(this->getPosition());
    Positioned::setPosition(position);
    Vector<3> newPos(this->getPosition());
    this->translateVertices(newPos - oldPos);
}

template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::translateVertices(const Vector<3> &translation) {
    auto *thisSpecific = static_cast<SpecificSolid*>(this);
    std::transform(thisSpecific->vertices.begin(), thisSpecific->vertices.end(), thisSpecific->vertices.begin(),
                   [&translation](const Vector<3> &vertex) {
                       return vertex + translation;
                   });
}

template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::calculateVerticesAndAxes() {
    auto *thisSpecific = static_cast<SpecificSolid*>(this);
    thisSpecific->vertices = this->applyOrientation(SpecificSolid::orientedVertices);
    thisSpecific->edgeAxes = this->applyOrientation(SpecificSolid::orientedEdgeAxes);
    thisSpecific->faceAxes = this->applyOrientation(SpecificSolid::orientedFaceAxes);
    Vector<3> pos(this->getPosition());
    this->translateVertices(pos);
}

template<typename SpecificSolid>
std::vector<std::string> PlatonicSolid<SpecificSolid>::getSupportedStrategies() const {
    return std::vector<std::string>{"sat", "tri-tri"};
}

template<typename SpecificSolid>
OverlapStrategy<3, 0> *PlatonicSolid<SpecificSolid>::createStrategy(const std::string &name) const {
    if (name == "sat")
        return new SATOverlap<SpecificSolid>;
    else if (name == "tri-tri")
        return new TriTriOverlap<SpecificSolid>;
    else
        return nullptr;
}

template<typename SpecificSolid>
std::string PlatonicSolid<SpecificSolid>::toWolfram() const {
    auto &thisSpecific = static_cast<const SpecificSolid&>(*this);
    auto triangles = thisSpecific.getTriangles();

    std::stringstream result;
    result << std::fixed;
    result << "Polygon[{";
    for (auto it = triangles.begin(); it != triangles.end() - 1; it++)
        result << "{" << (*it)[0] << ", " << (*it)[1] << ", " << (*it)[2] << "}," << std::endl << "    ";
    auto lastTri = triangles.back();
    result << "{" << lastTri[0] << ", " << lastTri[1] << ", " << lastTri[2] << "}}]";
    return result.str();
}
