//
// Created by PKua on 21.04.18.
//

#include "RegularOctahedron.h"

std::array<Vector<3>, 8> RegularOctahedron::orientedVertices;
std::array<Vector<3>, 4> RegularOctahedron::orientedFaceAxes;
std::array<Vector<3>, 6> RegularOctahedron::orientedEdgeAxes;

void RegularOctahedron::calculateStatic(const std::string &attr) {

}

double RegularOctahedron::projectionHalfsize(const Vector<3> &axis) const {
    return 0;
}

std::array<Vector<3>, 4> RegularOctahedron::getFaceAxes() const {
    return faceAxes;
}

std::array<Vector<3>, 6> RegularOctahedron::getEdgeAxes() const {
    return edgeAxes;
}

std::array<Vector<3>, 8> RegularOctahedron::getVertices() const {
    return this->vertices;
}

RegularOctahedron::RegularOctahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularOctahedron>{orientation} {
    this->calculateVertices();
    this->calculateAxes();
}
