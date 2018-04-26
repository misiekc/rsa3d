//
// Created by PKua on 21.04.18.
//

#include "RegularOctahedron.h"

std::array<Vector<3>, 4> RegularOctahedron::orientedFaceAxes;
std::array<Vector<3>, 6> RegularOctahedron::orientedEdgeAxes;

void RegularOctahedron::calculateStatic(const std::string &attr) {

}

double RegularOctahedron::projectionHalfsize(const Vector<3> &axis) const {
    return 0;
}
