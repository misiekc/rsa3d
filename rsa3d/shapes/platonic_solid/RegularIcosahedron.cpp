//
// Created by PKua on 21.04.18.
//

#include "RegularIcosahedron.h"

std::array<Vector<3>, 20> RegularIcosahedron::orientedVertices;
std::array<Vector<3>, 10> RegularIcosahedron::orientedFaceAxes;
std::array<Vector<3>, 15> RegularIcosahedron::orientedEdgeAxes;

std::array<Vector<3>, 20> RegularIcosahedron::getVertices() const {
    return this->vertices;
}

std::array<Vector<3>, 10> RegularIcosahedron::getFaceAxes() const {
    return this->faceAxes;
}

std::array<Vector<3>, 15> RegularIcosahedron::getEdgeAxes() const {
    return this->edgeAxes;
}

RegularIcosahedron::RegularIcosahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularIcosahedron>{orientation} {
    this->calculateVerticesAndAxes();
}

void RegularIcosahedron::calculateStatic(const std::string &attr) {

}

double RegularIcosahedron::projectionHalfsize(const Vector<3> &axis) const {
    return 0;
}
