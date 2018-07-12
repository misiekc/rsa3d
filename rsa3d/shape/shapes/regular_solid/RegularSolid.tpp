//
// Created by PKua on 21.04.18.
//

#include "../../../Vector.h"
#include "SATOverlap.h"
#include "TriTriOverlap.h"
#include <sstream>
#include <functional>
#include <iterator>
#include <fstream>

template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::orientedVertices;
template<typename SpecificSolid>
std::vector<std::vector<std::size_t>> RegularSolid<SpecificSolid>::orientedFaces;
template<typename SpecificSolid>
std::vector<std::array<std::size_t, 3>> RegularSolid<SpecificSolid>::orientedTriangles;
template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::orientedFaceAxes;
template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::orientedEdgeAxes;
template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::orientedVertexAxes;
template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::orientedMidedgeAxes;
template<typename SpecificSolid>
double RegularSolid<SpecificSolid>::normalizeFactor;
template<typename SpecificSolid>
double RegularSolid<SpecificSolid>::circumsphereRadius;
template<typename SpecificSolid>
double RegularSolid<SpecificSolid>::insphereRadius;


template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::initClass(const std::string &attr) {
    SpecificSolid::calculateStatic(attr);

    // No faces provided - generate wolfram notebook to recognize faces manually
    if (orientedFaces.empty()) {
        printNotebookWithVertices();
        exit(EXIT_SUCCESS);
    }

    discoverTriangles();
    normalizeVolume();
    calculateRadia();
    discoverAxes();

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
bool RegularSolid<SpecificSolid>::overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const {
    SpecificSolid other = dynamic_cast<const SpecificSolid&>(*s);     // Make a copy
    this->applyBC(bc, &other);
    SATOverlap<SpecificSolid> satOverlap;
    return satOverlap.overlap(this, &other);
}

template<typename SpecificSolid>
bool RegularSolid<SpecificSolid>::pointInside(BoundaryConditions<3> *bc, const Vector<3> &position,
                                              const Orientation<0> &orientation, double orientationRange) const {
    return (position - this->getPosition()).norm2() <= 4 * std::pow(SpecificSolid::insphereRadius, 2);
}

template<typename SpecificSolid>
bool RegularSolid<SpecificSolid>::isSeparatingAxis(const Vector<3> &axis, const SpecificSolid &other,
                                                    const Vector<3> &distance) const {
    auto &thisSpecific = static_cast<const SpecificSolid&>(*this);

    double distanceProj = std::abs(distance * axis);
    return distanceProj > thisSpecific.projectionHalfsize(axis) + other.projectionHalfsize(axis);
}

template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::store(std::ostream &f) const {
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
void RegularSolid<SpecificSolid>::restore(std::istream &f) {
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
const Matrix<3, 3> &RegularSolid<SpecificSolid>::getOrientationMatrix() const { return orientation; }

template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::setOrientationMatrix(const Matrix<3, 3> &orientation) {
    this->orientation = orientation;
}

template<typename SpecificSolid>
std::vector<std::string> RegularSolid<SpecificSolid>::getSupportedStrategies() const {
    return std::vector<std::string>{"sat", "tri-tri"};
}

template<typename SpecificSolid>
OverlapStrategy<3, 0> *RegularSolid<SpecificSolid>::createStrategy(const std::string &name) const {
    if (name == "sat")
        return new SATOverlap<SpecificSolid>;
    else if (name == "tri-tri")
        return new TriTriOverlap<SpecificSolid>;
    else
        return nullptr;
}

template<typename SpecificSolid>
std::string RegularSolid<SpecificSolid>::toWolfram() const {
    auto faces = this->getFaces();

    std::stringstream result;
    result << std::fixed;
    result << "Polygon[{ ";

    for (auto face = faces.begin(); face != faces.end(); face++) {
        result << "{";
        std::copy(face->begin(), face->end(), std::ostream_iterator<Vector<3>>(result, ", "));
        result.seekp(-2, std::ios_base::end);
        result << "}, " << std::endl;
    }
    result.seekp(-3, std::ios_base::end);
    result << " }]";

    return result.str();
}

template<typename SpecificSolid>
Shape<3, 0> *RegularSolid<SpecificSolid>::clone() const {
    auto &thisSpecific = static_cast<const SpecificSolid&>(*this);
    return new SpecificSolid(thisSpecific);
}

template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::applyOrientation(const std::vector<Vector<3>> &vectors) const {
    std::vector<Vector<3>> result;
    result.reserve(vectors.size());
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(result), [this](const Vector<3> &vector) {
        return this->orientation * vector;
    });
    return result;
}

template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::applyPosition(const std::vector<Vector<3>> &vectors) const {
    std::vector<Vector<3>> result;
    result.reserve(vectors.size());
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(result), [this](const Vector<3> &vector) {
        return this->orientation * vector + this->getPosition();
    });
    return result;
}

template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::getVertices() const { return applyPosition(orientedVertices); }

template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::getEdgeAxes() const { return applyOrientation(orientedEdgeAxes); }

template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::getFaceAxes() const { return applyOrientation(orientedFaceAxes); }

template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::getVertexAxes() const { return applyOrientation(orientedVertexAxes); }

template<typename SpecificSolid>
std::vector<Vector<3>> RegularSolid<SpecificSolid>::getMidegdeAxes() const { return applyOrientation(orientedMidedgeAxes); }

template<typename SpecificSolid>
intersection::tri_polyh RegularSolid<SpecificSolid>::getTriangles() const {
    intersection::tri_polyh result;
    result.reserve(orientedTriangles.size());

    auto vertices = this->getVertices();
    std::transform(orientedTriangles.begin(), orientedTriangles.end(), std::back_inserter(result),
                   [&vertices](const std::array<std::size_t, 3> &tri) {
                       return intersection::tri3D{{vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]}};
                   });
    return result;
}

template<typename SpecificSolid>
intersection::face_polyh RegularSolid<SpecificSolid>::getFaces() const {
    intersection::face_polyh result;
    result.reserve(orientedFaces.size());

    auto vertices = this->getVertices();
    std::transform(orientedFaces.begin(), orientedFaces.end(), std::back_inserter(result),
                   [&vertices](const std::vector<std::size_t> &face) {
                       intersection::face3D faceResult;
                       faceResult.reserve(face.size());
                       std::transform(face.begin(), face.end(), std::back_inserter(faceResult), [&vertices](size_t i) {
                           return vertices[i];
                       });
                       return faceResult;
                   });
    return result;
}

template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::discoverTriangles() {
    for (const auto &face : orientedFaces) {
        std::size_t numVertices = face.size();
        if (numVertices < 3)
            throw std::runtime_error("number of vertices in face < 3");
        for (std::size_t i = 2; i < numVertices; i++)
            orientedTriangles.push_back({{face[0], face[i - 1], face[i]}});
    }
    std::cout << "[RegularSolid::discoverTriangles] Faces: " << orientedFaces.size();
    std::cout << ", triangles: " << orientedTriangles.size() << std::endl;
}

template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::normalizeVolume() {
    double volume = std::accumulate(orientedTriangles.begin(), orientedTriangles.end(), 0.,
            [](double volume, const std::array<std::size_t, 3> &triIdx) {
        intersection::tri3D tri{{orientedVertices[triIdx[0]], orientedVertices[triIdx[1]], orientedVertices[triIdx[2]]}};
        double tetrahedronVolume = 1./6 * tri[0] * ((tri[1] - tri[0]) ^ (tri[2] - tri[0]));
        return volume + std::abs(tetrahedronVolume);
    });
    normalizeFactor = std::pow(volume, -1./3);
    std::transform(orientedVertices.begin(), orientedVertices.end(), orientedVertices.begin(),
                   [](const Vector<3> &vertex) {
        return vertex * normalizeFactor;
    });
    std::cout << "[RegularSolid::normalizeVolume] Normalizing factor: " << normalizeFactor << std::endl;
}

template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::calculateRadia() {
    circumsphereRadius = 0;
    for (const auto &vertex : orientedVertices) {
        double vertexDistance = vertex.norm();
        if (vertexDistance > circumsphereRadius)
            circumsphereRadius = vertexDistance;
    }

    insphereRadius = std::numeric_limits<double>::infinity();
    for (const auto &face : orientedFaces) {
        Vector<3> edge1 = orientedVertices[face[1]] - orientedVertices[face[0]];
        Vector<3> edge2 = orientedVertices[face[2]] - orientedVertices[face[0]];
        Vector<3> normalAxis = edge1 ^ edge2;
        normalAxis /= normalAxis.norm();
        double faceDistance = std::abs(orientedVertices[face[0]] * normalAxis);
        if (faceDistance < insphereRadius)
            insphereRadius = faceDistance;
    }

    std::cout << "[RegularSolid::calculateRadia] Circumsphere radius: " << circumsphereRadius;
    std::cout << ", insphere radius: " << insphereRadius << std::endl;
}

template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::discoverAxes() {
    for (const auto &vertex : orientedVertices) {
        Vector<3> axis = vertex / vertex.norm();
        addUniqueAxis(orientedVertexAxes, axis);
    }

    for (const auto &faceIdx : orientedFaces) {
        intersection::face3D face(faceIdx.size());
        std::transform(faceIdx.begin(), faceIdx.end(), face.begin(), [](std::size_t idx) {
            return orientedVertices[idx];
        });

        Vector<3> faceAxis = (face[1] - face[0]) ^ (face[2] - face[0]);
        faceAxis /= faceAxis.norm();
        addUniqueAxis(orientedFaceAxes, faceAxis);

        for (std::size_t i = 0; i < face.size(); i++) {
            Vector<3> edgeAxis = face[(i + 1) % face.size()] - face[i];
            edgeAxis = edgeAxis / edgeAxis.norm();
            addUniqueAxis(orientedEdgeAxes, edgeAxis);

            Vector<3> midedgeAxis = (face[(i + 1) % face.size()] + face[i]) / 2.;
            midedgeAxis = midedgeAxis / midedgeAxis.norm();
            addUniqueAxis(orientedMidedgeAxes, midedgeAxis);
        }
    }

    std::cout << "[RegularSolid::discoverAxes] Discovered ";
    std::cout << orientedFaceAxes.size() << " face axes, ";
    std::cout << orientedEdgeAxes.size() << " edge axes, ";
    std::cout << orientedVertexAxes.size() << " vertex axes, ";
    std::cout << orientedMidedgeAxes.size() << " midedge axes" << std::endl;
}

template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::addUniqueAxis(std::vector<Vector<3>> &axes, const Vector<3> &newAxis) {
    auto axisCompare = [&newAxis](const Vector<3> axis) {
        double product = std::abs(axis * newAxis);
        return std::abs(product - 1) < 0.000000001;
    };
    auto it = find_if(axes.begin(), axes.end(), axisCompare);
    if (it == axes.end())
        axes.push_back(newAxis);
}

template<typename SpecificSolid>
void RegularSolid<SpecificSolid>::printNotebookWithVertices() {
    std::ofstream file("goodluck.nb");
    if (!file)  die("Cannot open goodluck.nb to store vertices");

    file << "vertices={";
    std::copy(orientedVertices.begin(), orientedVertices.end(), std::ostream_iterator<const Vector<3>>(file, ", "));
    file.seekp(-2, std::ios_base::end);
    file << "};" << std::endl;
    file << "labels=Table[Text[Style[i-1,Large],vertices[[i]]],{i,1," << orientedVertices.size() << "}];" << std::endl << std::endl;
    file << "(* List here faces in a format of {{vertexIdx1, vertexIdx2, ...}, {face2}, {face3}, ...} *)" << std::endl;
    file << "faces={};" << std::endl << std::endl;
    file << "facesCoord=Map[vertices[[#+1]]&,faces,{2}];" << std::endl;
    file << "Graphics3D[Join[labels,{Polygon[facesCoord]}]]";

    file.close();

    std::cout << "[RegularSolid::printNotebookWithVertices] No faces provided. Vertices printed to goodluck.nb. ";
    std::cout << "Use it to recognize faces manually." << std::endl;
}