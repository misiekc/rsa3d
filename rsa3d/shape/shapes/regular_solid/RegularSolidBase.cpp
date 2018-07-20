//
// Created by PKua on 20.07.18.
//

#include <iterator>
#include <fstream>
#include "RegularSolidBase.h"
#include "SATOverlapRS.h"
#include "TriTriOverlapRS.h"


std::vector<Vector<3>> RegularSolidBase::orientedVertices;
std::vector<std::vector<std::size_t>> RegularSolidBase::orientedFaces;
std::vector<std::array<std::size_t, 3>> RegularSolidBase::orientedTriangles;
std::vector<Vector<3>> RegularSolidBase::orientedFaceAxes;
std::vector<Vector<3>> RegularSolidBase::orientedEdgeAxes;
std::vector<Vector<3>> RegularSolidBase::orientedVertexAxes;
std::vector<Vector<3>> RegularSolidBase::orientedMidedgeAxes;
double RegularSolidBase::normalizeFactor;
double RegularSolidBase::circumsphereRadius;
double RegularSolidBase::insphereRadius;

void RegularSolidBase::initClass(const std::string &attr) {
    // No faces provided - generate wolfram notebook to recognize faces manually
    if (orientedFaces.empty()) {
        printFaceHelperNotebook();
        exit(EXIT_SUCCESS);
    }

    discoverTriangles();
    normalizeVolume();
    calculateRadia();
    discoverAxes();

    Shape::setNeighbourListCellSize(2*circumsphereRadius);
    Shape::setVoxelSpatialSize(2*insphereRadius/std::sqrt(3.));
}

bool RegularSolidBase::pointInside(BoundaryConditions<3> *bc, const Vector<3> &position,
                                              const Orientation<0> &orientation, double orientationRange) const {
    Vector<3> thisPos = this->getPosition();
    Vector<3> pointDelta = position - thisPos + bc->getTranslation(thisPos, position);
    return pointDelta.norm2() <= 4*insphereRadius*insphereRadius;
}

void RegularSolidBase::store(std::ostream &f) const {
    Shape::store(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            d = this->orientation(i, j);
            f.write((char *)(&d), sizeof(double));
        }
    }
}

void RegularSolidBase::restore(std::istream &f) {
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

const Matrix<3, 3> &RegularSolidBase::getOrientationMatrix() const { return orientation; }

void RegularSolidBase::setOrientationMatrix(const Matrix<3, 3> &orientation) {
    this->orientation = orientation;
}

std::vector<std::string> RegularSolidBase::getSupportedStrategies() const {
    return std::vector<std::string>{"sat", "tri-tri"};
}

OverlapStrategy<3, 0> *RegularSolidBase::createStrategy(const std::string &name) const {
    if (name == "sat")
        return new SATOverlapRS{};
    else if (name == "tri-tri")
        return new TriTriOverlapRS{};
    else
        return nullptr;
}

std::string RegularSolidBase::toWolfram() const {
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

std::vector<Vector<3>> RegularSolidBase::applyOrientation(const std::vector<Vector<3>> &vectors) const {
    std::vector<Vector<3>> result;
    result.reserve(vectors.size());
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(result), [this](const Vector<3> &vector) {
        return this->orientation * vector;
    });
    return result;
}

std::vector<Vector<3>> RegularSolidBase::applyPosition(const std::vector<Vector<3>> &vectors) const {
    std::vector<Vector<3>> result;
    result.reserve(vectors.size());
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(result), [this](const Vector<3> &vector) {
        return this->orientation * vector + this->getPosition();
    });
    return result;
}

std::vector<Vector<3>> RegularSolidBase::getVertices() const { return applyPosition(orientedVertices); }
std::vector<Vector<3>> RegularSolidBase::getEdgeAxes() const { return applyOrientation(orientedEdgeAxes); }
std::vector<Vector<3>> RegularSolidBase::getFaceAxes() const { return applyOrientation(orientedFaceAxes); }
std::vector<Vector<3>> RegularSolidBase::getVertexAxes() const { return applyOrientation(orientedVertexAxes); }
std::vector<Vector<3>> RegularSolidBase::getMidegdeAxes() const { return applyOrientation(orientedMidedgeAxes); }

intersection::tri_polyh RegularSolidBase::getTriangles() const {
    auto vertices = this->getVertices();
    auto verticesIdxToTri = [&vertices](const std::array<std::size_t, 3> &tri) {
        return intersection::tri3D{{vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]}};
    };

    intersection::tri_polyh triangles;
    triangles.reserve(orientedTriangles.size());
    std::transform(orientedTriangles.begin(), orientedTriangles.end(), std::back_inserter(triangles), verticesIdxToTri);
    return triangles;
}

intersection::face_polyh RegularSolidBase::getFaces() const {
    auto vertices = this->getVertices();
    auto verticesIdxToFace = [&vertices](const std::vector<std::size_t> &face) {
        intersection::face3D faceResult;
        faceResult.reserve(face.size());
        auto idxToVertex = [&vertices](size_t i) { return vertices[i]; };
        std::transform(face.begin(), face.end(), std::back_inserter(faceResult), idxToVertex);
        return faceResult;
    };

    intersection::face_polyh faces;
    faces.reserve(orientedFaces.size());
    std::transform(orientedFaces.begin(), orientedFaces.end(), std::back_inserter(faces), verticesIdxToFace);
    return faces;
}

void RegularSolidBase::discoverTriangles() {
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

void RegularSolidBase::normalizeVolume() {
    auto triangleVolumeAcumulator = [](double volume, const std::array<std::size_t, 3> &triIdx) {
        intersection::tri3D tri{{orientedVertices[triIdx[0]], orientedVertices[triIdx[1]], orientedVertices[triIdx[2]]}};
        double tetrahedronVolume = 1./6 * tri[0] * ((tri[1] - tri[0]) ^ (tri[2] - tri[0]));
        return volume + std::abs(tetrahedronVolume);
    };
    double volume = std::accumulate(orientedTriangles.begin(), orientedTriangles.end(), 0., triangleVolumeAcumulator);

    normalizeFactor = std::pow(volume, -1./3);
    auto normalizeVertex = [](const Vector<3> &vertex) { return vertex * normalizeFactor; };
    std::transform(orientedVertices.begin(), orientedVertices.end(), orientedVertices.begin(), normalizeVertex);
    std::cout << "[RegularSolid::normalizeVolume] Normalizing factor: " << normalizeFactor << std::endl;
}

void RegularSolidBase::calculateRadia() {
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

void RegularSolidBase::discoverAxes() {
    for (const auto &vertex : orientedVertices) {
        Vector<3> axis = vertex / vertex.norm();
        addUniqueAxis(orientedVertexAxes, axis);
    }

    for (const auto &faceIdx : orientedFaces) {
        intersection::face3D face(faceIdx.size());
        auto idxToVertex = [](std::size_t idx) { return orientedVertices[idx]; };
        std::transform(faceIdx.begin(), faceIdx.end(), face.begin(), idxToVertex);

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

void RegularSolidBase::addUniqueAxis(std::vector<Vector<3>> &axes, const Vector<3> &newAxis) {
    auto isNewAxisUnique = [&newAxis](const Vector<3> axis) {
        double product = std::abs(axis * newAxis);
        return std::abs(product - 1) < 0.000000001;
    };
    auto it = find_if(axes.begin(), axes.end(), isNewAxisUnique);
    if (it == axes.end())
        axes.push_back(newAxis);
}

void RegularSolidBase::printFaceHelperNotebook() {
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

    std::cout << "[RegularSolid::printFaceHelperNotebook] No faces provided. Vertices printed to goodluck.nb. ";
    std::cout << "Use it to recognize faces manually." << std::endl;
}

double RegularSolidBase::projectionHalfsize(const Vector<3> &axis) const {
    throw std::runtime_error("unimplemented");
}
