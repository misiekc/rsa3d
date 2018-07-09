//
// Created by PKua on 21.04.18.
//

#include "RegularOctahedron.h"

void RegularOctahedron::calculateStatic(const std::string &attr) {
    auto &vertices = PlatonicSolid<RegularOctahedron>::orientedVertices;
    auto &edgeAxes = PlatonicSolid<RegularOctahedron>::orientedEdgeAxes;
    auto &faceAxes = PlatonicSolid<RegularOctahedron>::orientedFaceAxes;
    
    vertices = {Vector<3>{{ 1,  0,  0}} * edgeFactor,
                         Vector<3>{{-1,  0,  0}} * edgeFactor,
                         Vector<3>{{ 0,  1,  0}} * edgeFactor,
                         Vector<3>{{ 0, -1,  0}} * edgeFactor,
                         Vector<3>{{ 0,  0,  1}} * edgeFactor,
                         Vector<3>{{ 0,  0, -1}} * edgeFactor};

    edgeAxes = {vertices[4] - vertices[2],
                         vertices[4] - vertices[1],
                         vertices[3] - vertices[4],
                         vertices[0] - vertices[4],
                         vertices[0] - vertices[3],
                         vertices[2] - vertices[0]};

    faceAxes = {edgeAxes[3] ^ edgeAxes[0],
                         edgeAxes[0] ^ edgeAxes[1],
                         edgeAxes[2] ^ edgeAxes[1],
                         edgeAxes[2] ^ edgeAxes[3]};
}

double RegularOctahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xHalfsize = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yHalfsize = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zHalfsize = std::abs(this->getOrientationMatrix().column(2) * axis);

    return std::max(std::max(xHalfsize, yHalfsize), zHalfsize) * edgeFactor;
}

RegularOctahedron::RegularOctahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularOctahedron>{orientation} {
    this->calculateVerticesAndAxes();
}

intersection::tri_polyh RegularOctahedron::getTriangles() const {
    auto vertices = this->getVertices();
    return intersection::tri_polyh{
            {{vertices[0], vertices[2], vertices[4]}},
            {{vertices[1], vertices[5], vertices[3]}},
            {{vertices[4], vertices[2], vertices[1]}},
            {{vertices[3], vertices[5], vertices[0]}},
            {{vertices[4], vertices[1], vertices[3]}},
            {{vertices[2], vertices[0], vertices[5]}},
            {{vertices[3], vertices[0], vertices[4]}},
            {{vertices[1], vertices[2], vertices[5]}}
    };
}
