//
// Created by PKua on 20.07.18.
//

#include <iterator>
#include <fstream>
#include <iomanip>
#include "RegularSolidBase.h"
#include "SATOverlapRS.h"
#include "TriTriOverlapRS.h"


//#define PRINT_FACE_AXES
//#define PRINT_VERTEX_AXES


void RegularSolidBase::complementShapeData(ShapeData &shapeDataToComplement) {
    // No faces provided - generate wolfram notebook to recognize faces manually
    if (shapeDataToComplement.orientedFaces.empty()) {
        printFaceHelperNotebook(shapeDataToComplement);
        exit(EXIT_SUCCESS);
    }

    double circumsphereRadius, insphereRadius;

    normalizeFacesOrientation(shapeDataToComplement);
    discoverEdges(shapeDataToComplement);
    discoverTriangles(shapeDataToComplement);
    normalizeVolume(shapeDataToComplement);
    calculateRadia(shapeDataToComplement, circumsphereRadius, insphereRadius);
    discoverAxes(shapeDataToComplement);

    ShapeStaticInfo<3, 0> shapeInfo;
    shapeInfo.setCircumsphereRadius(circumsphereRadius);
    shapeInfo.setInsphereRadius(insphereRadius);
    shapeInfo.setExclusionZoneMaxSpan(circumsphereRadius + insphereRadius);
    shapeInfo.setCreateShapeImpl([](auto dummy) -> Shape* { return nullptr; });   // To be replaced in RegularSolid.tpp
    Shape::setShapeStaticInfo(shapeInfo);

#ifndef NDEBUG
    shapeDataToComplement.printInfo(std::cout);
#endif
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

    for (auto &face : faces) {
        result << "{";
        std::copy(face.begin(), face.end(), std::ostream_iterator<Vector<3>>(result, ", "));
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
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(result), [this](auto &vector) {
        return this->orientation * vector;
    });
    return result;
}

std::vector<Vector<3>> RegularSolidBase::applyPosition(const std::vector<Vector<3>> &vectors) const {
    std::vector<Vector<3>> result;
    result.reserve(vectors.size());
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(result), [this](auto &vector) {
        return this->orientation * vector + this->getPosition();
    });
    return result;
}

bool RegularSolidBase::pointInside(BoundaryConditions<3> *bc, const Vector<3> &position,
                                   const Orientation<0> &orientation, double orientationRange) const {
    Vector<3> bcPos = position + bc->getTranslation(this->getPosition(), position);
    auto vertices = this->getVertices();
    switch (pointInsideFace(bcPos, vertices)) {
        case TRUE:      return true;
        case FALSE:     return false;
        default:        return pointInsideEdge(bcPos, vertices) || pointInsideVertex(bcPos, vertices);
    }
}

inline RegularSolidBase::PIResult
RegularSolidBase::pointInsideFace(const Vector<3> &point, const std::vector<Vector<3>> &vertices) const {
    auto faceNormals = this->getFaceNormals();
    bool insideSolid = true;
    for (std::size_t faceI = 0; faceI < this->shapeData->orientedFaces.size(); faceI++) {
        const auto &face = this->shapeData->orientedFaces[faceI];
        Vector<3> faceNormal = faceNormals[faceI];
        Vector<3> firstVertex = vertices[face[0]];
        double faceDist = (point - firstVertex) * faceNormal;
        if (faceDist > 0) {
            insideSolid = false;
            if (faceDist > Shape::getInsphereRadius())
                return FALSE;
            else if (projectionInsideFace(point, vertices, face, faceNormal))
                return TRUE;
        }
    }

    if (insideSolid)
        return TRUE;

    return UNKNOWN;
}

inline bool RegularSolidBase::projectionInsideFace(const Vector<3> &point, const std::vector<Vector<3>> &vertices,
                                                   const std::vector<size_t> &face, const Vector<3> &faceNormal) const {
    for (std::size_t vertexI = 0; vertexI < face.size(); vertexI++) {
        Vector<3> edgeBeg = vertices[face[vertexI]];
        Vector<3> edgeEnd = vertices[face[(vertexI + 1) % face.size()]];
        Vector<3> edgeVec = edgeEnd - edgeBeg;
        Vector<3> posVec = point - edgeBeg;
        if ((edgeVec ^ faceNormal) * posVec > 0)
            return false;
    }
    return true;
}

inline bool RegularSolidBase::pointInsideEdge(const Vector<3> &point, const std::vector<Vector<3>> &vertices) const {
    for (const auto &edge : this->shapeData->orientedEdges) {
        Vector<3> edgeBeg = vertices[edge.first()];
        Vector<3> edgeEnd = vertices[edge.second()];
        Vector<3> edgeVec = edgeEnd - edgeBeg;
        Vector<3> posFromBeg = point - edgeBeg;

        Vector<3> orthPosComponent = posFromBeg - posFromBeg.projectOn(edgeVec);
        if (orthPosComponent.norm2() <= std::pow(Shape::getInsphereRadius(), 2)) {
            Vector<3> posFromEnd = point - edgeEnd;
            if (posFromBeg * edgeVec >= 0 && posFromEnd * edgeVec <= 0)
                return true;
        }
    }
    return false;
}

inline bool RegularSolidBase::pointInsideVertex(const Vector<3> &point, const std::vector<Vector<3>> &vertices) const {
    for (const auto &vertex : vertices)
        if ((point - vertex).norm2() <= std::pow(Shape::getInsphereRadius(), 2))
            return true;
    return false;
}

std::vector<Vector<3>> RegularSolidBase::getVertices() const {
    return this->applyPosition(this->shapeData->orientedVertices);
}

std::vector<Vector<3>> RegularSolidBase::getEdgeAxes() const {
    return this->applyOrientation(this->shapeData->orientedEdgeAxes);
}

std::vector<Vector<3>> RegularSolidBase::getFaceAxes() const {
    return this->applyOrientation(this->shapeData->orientedFaceAxes);
}

std::vector<Vector<3>> RegularSolidBase::getFaceNormals() const {
    return this->applyOrientation(this->shapeData->orientedFaceNormals);
}

std::vector<Vector<3>> RegularSolidBase::getVertexAxes() const {
    return this->applyOrientation(this->shapeData->orientedVertexAxes);
}

std::vector<Vector<3>> RegularSolidBase::getMidegdeAxes() const {
    return this->applyOrientation(this->shapeData->orientedMidedgeAxes);
}

PolyhedronTriangulation RegularSolidBase::getTriangulation() const {
    auto vertices = this->getVertices();
    auto verticesIdxToTri = [&vertices](const std::array<std::size_t, 3> &tri) {
        return Triangle3D{{vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]}};
    };

    PolyhedronTriangulation triangles;
    triangles.reserve(this->shapeData->orientedTriangles.size());
    std::transform(this->shapeData->orientedTriangles.begin(), this->shapeData->orientedTriangles.end(),
                   std::back_inserter(triangles), verticesIdxToTri);
    return triangles;
}

PolyhedronFaces RegularSolidBase::getFaces() const {
    auto vertices = this->getVertices();
    auto verticesIdxToFace = [&vertices](const std::vector<std::size_t> &face) {
        Face3D faceResult;
        faceResult.reserve(face.size());
        auto idxToVertex = [&vertices](size_t i) { return vertices[i]; };
        std::transform(face.begin(), face.end(), std::back_inserter(faceResult), idxToVertex);
        return faceResult;
    };

    PolyhedronFaces faces;
    faces.reserve(this->shapeData->orientedFaces.size());
    std::transform(this->shapeData->orientedFaces.begin(), this->shapeData->orientedFaces.end(),
                   std::back_inserter(faces), verticesIdxToFace);
    return faces;
}

Face3D RegularSolidBase::ShapeData::faceFromVertexIdx(const std::vector<size_t> &vertexIdx) const {
    Face3D face(vertexIdx.size());
    auto idxToVertex = [this](std::size_t idx) { return this->orientedVertices[idx]; };
    transform(vertexIdx.begin(), vertexIdx.end(), face.begin(), idxToVertex);
    return face;
}

void RegularSolidBase::normalizeFacesOrientation(ShapeData &shapeDataToComplement) {
    for (auto &faceVertexIdx : shapeDataToComplement.orientedFaces) {
        Face3D face = shapeDataToComplement.faceFromVertexIdx(faceVertexIdx);
        Vector<3> faceAxis = (face[1] - face[0]) ^ (face[2] - face[0]);
        if (faceAxis * face[0] < 0)
            std::reverse(faceVertexIdx.begin(), faceVertexIdx.end());
    }
}

void RegularSolidBase::discoverEdges(ShapeData &shapeDataToComplement) {
    for (const auto &face : shapeDataToComplement.orientedFaces) {
        for (std::size_t vertexI = 0; vertexI < face.size(); vertexI++) {
            Edge edge(face[(vertexI + 1) % face.size()], face[vertexI]);
            if (std::find(shapeDataToComplement.orientedEdges.begin(), shapeDataToComplement.orientedEdges.end(), edge)
                == shapeDataToComplement.orientedEdges.end())
            {
                shapeDataToComplement.orientedEdges.push_back(edge);
            }
        }
    }
}

void RegularSolidBase::discoverTriangles(ShapeData &shapeDataToComplement) {
    for (const auto &face : shapeDataToComplement.orientedFaces) {
        std::size_t numVertices = face.size();
        if (numVertices < 3)
            throw std::runtime_error("number of vertices in face < 3");
        for (std::size_t i = 2; i < numVertices; i++)
            shapeDataToComplement.orientedTriangles.push_back({{face[0], face[i - 1], face[i]}});
    }
}

void RegularSolidBase::normalizeVolume(ShapeData &shapeDataToComplement) {
    auto triangleVolumeAcumulator = [shapeDataToComplement](double volume, const std::array<std::size_t, 3> &triIdx) {
        Triangle3D tri{{shapeDataToComplement.orientedVertices[triIdx[0]],
                        shapeDataToComplement.orientedVertices[triIdx[1]],
                        shapeDataToComplement.orientedVertices[triIdx[2]]}};
        double tetrahedronVolume = 1./6 * tri[0] * ((tri[1] - tri[0]) ^ (tri[2] - tri[0]));
        return volume + std::abs(tetrahedronVolume);
    };
    double volume = std::accumulate(shapeDataToComplement.orientedTriangles.begin(),
                                    shapeDataToComplement.orientedTriangles.end(), 0., triangleVolumeAcumulator);

    shapeDataToComplement.normalizeFactor = std::pow(volume, -1./3);
    auto normalizeVertex = [shapeDataToComplement](const Vector<3> &vertex) {
        return vertex * shapeDataToComplement.normalizeFactor;
    };
    std::transform(shapeDataToComplement.orientedVertices.begin(), shapeDataToComplement.orientedVertices.end(),
                   shapeDataToComplement.orientedVertices.begin(), normalizeVertex);
}

void RegularSolidBase::calculateRadia(const ShapeData &shapeData, double &circumsphereRadius,
                                      double &insphereRadius) {
    circumsphereRadius = 0;
    for (const auto &vertex : shapeData.orientedVertices) {
        double vertexDistance = vertex.norm();
        if (vertexDistance > circumsphereRadius)
            circumsphereRadius = vertexDistance;
    }

    insphereRadius = std::numeric_limits<double>::infinity();
    for (const auto &face : shapeData.orientedFaces) {
        Vector<3> edge1 = shapeData.orientedVertices[face[1]]
                          - shapeData.orientedVertices[face[0]];
        Vector<3> edge2 = shapeData.orientedVertices[face[2]]
                          - shapeData.orientedVertices[face[0]];
        Vector<3> normalAxis = edge1 ^ edge2;
        normalAxis /= normalAxis.norm();
        double faceDistance = std::abs(shapeData.orientedVertices[face[0]] * normalAxis);
        if (faceDistance < insphereRadius)
            insphereRadius = faceDistance;
    }
}

void RegularSolidBase::discoverAxes(ShapeData &shapeDataToComplement) {
    for (const auto &vertex : shapeDataToComplement.orientedVertices)
        addUniqueAxis(shapeDataToComplement.orientedVertexAxes, vertex.normalized());

    for (const auto &faceVertexIdx : shapeDataToComplement.orientedFaces) {
        Face3D face = shapeDataToComplement.faceFromVertexIdx(faceVertexIdx);

        Vector<3> faceAxis = ((face[1] - face[0]) ^ (face[2] - face[0])).normalized();
        addUniqueAxis(shapeDataToComplement.orientedFaceAxes, faceAxis);
        shapeDataToComplement.orientedFaceNormals.push_back(faceAxis);

        for (std::size_t i = 0; i < face.size(); i++) {
            Vector<3> edgeAxis = face[(i + 1) % face.size()] - face[i];
            addUniqueAxis(shapeDataToComplement.orientedEdgeAxes, edgeAxis.normalized());

            Vector<3> midedgeAxis = (face[(i + 1) % face.size()] + face[i]) / 2.;
            addUniqueAxis(shapeDataToComplement.orientedMidedgeAxes, midedgeAxis.normalized());
        }
    }
}

void RegularSolidBase::ShapeData::printInfo(std::ostream &out) const {
    out << "[RegularSolid::initClass] Faces: " << orientedFaces.size();
    out << ", triangles: " << orientedTriangles.size() << std::endl;
    out << "[RegularSolid::initClass] Normalizing factor: " << normalizeFactor << std::endl;
    out << "[RegularSolid::initClass] Circumsphere radius: " << Shape::getCircumsphereRadius();
    out << ", insphere radius: " << Shape::getInsphereRadius() << std::endl;
    out << "[RegularSolid::initClass] Discovered ";
    out << orientedFaceAxes.size() << " face axes, ";
    out << orientedEdgeAxes.size() << " edge axes, ";
    out << orientedVertexAxes.size() << " vertex axes, ";
    out << orientedMidedgeAxes.size() << " midedge axes" << std::endl;

#ifdef PRINT_FACE_AXES
    out << "[RegularSolid::initClass] Face axes:" << std::endl;
    this->printAxes(orientedFaceAxes, out);
#endif

#ifdef PRINT_VERTEX_AXES
    out << "[RegularSolid::initClass] Vertex axes:" << std::endl;
    this->printAxes(orientedVertexAxes,out);
#endif
}

void RegularSolidBase::ShapeData::printAxes(const std::vector<Vector<3>> &axes, std::ostream &out) const {
    std::size_t index = 1;
    for (const auto &axis : axes)
        out << std::setw(3) << std::right  << (index++) << ".: " << axis << std::endl;
}

void RegularSolidBase::addUniqueAxis(std::vector<Vector<3>> &axes, const Vector<3> &newAxis) {
    auto isNewAxisUnique = [&newAxis](const Vector<3> &axis) {
        double product = std::abs(axis * newAxis);
        return std::abs(product - 1) < 0.000000001;
    };
    auto it = find_if(axes.begin(), axes.end(), isNewAxisUnique);
    if (it == axes.end())
        axes.push_back(newAxis);
}

void RegularSolidBase::printFaceHelperNotebook(const ShapeData &shapeData) {
    std::ofstream file("goodluck.nb");
    if (!file)  die("Cannot open goodluck.nb to store vertices");

    file << "vertices={";
    std::copy(shapeData.orientedVertices.begin(), shapeData.orientedVertices.end(),
              std::ostream_iterator<const Vector<3>>(file, ", "));
    file.seekp(-2, std::ios_base::end);
    file << "};" << std::endl;
    file << "labels=Table[Text[Style[i-1,Large],vertices[[i]]],{i,1," << shapeData.orientedVertices.size();
    file << "}];" << std::endl << std::endl;
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
