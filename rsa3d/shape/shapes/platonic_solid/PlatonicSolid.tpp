//
// Created by PKua on 21.04.18.
//

#include "../../../Vector.h"
#include "SATOverlap.h"
#include "TriTriOverlap.h"
#include <sstream>

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::orientedVertices;

template<typename SpecificSolid>
intersection::face_polyh PlatonicSolid<SpecificSolid>::orientedFaces;

template<typename SpecificSolid>
intersection::tri_polyh PlatonicSolid<SpecificSolid>::orientedTriangles;

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::orientedFaceAxes;

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::orientedEdgeAxes;

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::orientedVertexAxes;

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::orientedMidedgeAxes;


template<typename SpecificSolid>
void PlatonicSolid<SpecificSolid>::initClass(const std::string &attr) {
    SpecificSolid::calculateStatic(attr);
    normalizeAxes();

    Shape::setNeighbourListCellSize(2 * SpecificSolid::circumsphereRadius);
    Shape::setVoxelSpatialSize(2 * SpecificSolid::insphereRadius / std::sqrt(3.));

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
bool PlatonicSolid<SpecificSolid>::overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const {
    SpecificSolid other = dynamic_cast<const SpecificSolid&>(*s);     // Make a copy
    this->applyBC(bc, &other);
    SATOverlap<SpecificSolid> satOverlap;
    return satOverlap.overlap(this, &other);
}

template<typename SpecificSolid>
bool PlatonicSolid<SpecificSolid>::pointInside(BoundaryConditions<3> *bc, const Vector<3> &position,
                                              const Orientation<0> &orientation, double orientationRange) const {
    return (position - this->getPosition()).norm2() <= 4 * std::pow(SpecificSolid::insphereRadius, 2);
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
void PlatonicSolid<SpecificSolid>::calculateVerticesAndAxes() {

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

template<typename SpecificSolid>
Shape<3, 0> *PlatonicSolid<SpecificSolid>::clone() const {
    auto &thisSpecific = static_cast<const SpecificSolid&>(*this);
    return new SpecificSolid(thisSpecific);
}

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::applyOrientation(const std::vector<Vector<3>> &vectors) const {
    std::vector<Vector<3>> result;
    result.reserve(vectors.size());
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(result), [this](const Vector<3> &vector) {
        return this->orientation * vector;
    });
    return result;
}

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::applyPosition(const std::vector<Vector<3>> &vectors) const {
    std::vector<Vector<3>> result;
    result.reserve(vectors.size());
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(result), [this](const Vector<3> &vector) {
        return this->orientation * vector + this->getPosition();
    });
    return result;
}

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::getVertices() const { return applyPosition(orientedVertices); }

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::getEdgeAxes() const { return applyOrientation(orientedEdgeAxes); }

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::getFaceAxes() const { return applyOrientation(orientedFaceAxes); }

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::getVertexAxes() const { return applyOrientation(orientedVertexAxes); }

template<typename SpecificSolid>
std::vector<Vector<3>> PlatonicSolid<SpecificSolid>::getMidegdeAxes() const { return applyOrientation(orientedMidedgeAxes); }