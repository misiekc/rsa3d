//
// Created by PKua on 21.04.18.
//

#include "RegularDodecahedron.h"

std::array<Vector<3>, 12> RegularDodecahedron::orientedVertices;
std::array<Vector<3>, 6> RegularDodecahedron::orientedFaceAxes;
std::array<Vector<3>, 15> RegularDodecahedron::orientedEdgeAxes;

std::array<Vector<3>, 12> RegularDodecahedron::getVertices() const {
    return this->vertices;
}

std::array<Vector<3>, 6> RegularDodecahedron::getFaceAxes() const {
    return this->faceAxes;
}

std::array<Vector<3>, 15> RegularDodecahedron::getEdgeAxes() const {
    return this->edgeAxes;
}

RegularDodecahedron::RegularDodecahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularDodecahedron>{orientation} {
    this->calculateVerticesAndAxes();
}

void RegularDodecahedron::calculateStatic(const std::string &attr) {

}

double RegularDodecahedron::projectionHalfsize(const Vector<3> &axis) const {
    return 0;
}
