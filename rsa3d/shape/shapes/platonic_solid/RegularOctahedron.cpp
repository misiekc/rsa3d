//
// Created by PKua on 21.04.18.
//

#include "RegularOctahedron.h"

std::array<Vector<3>, 8> RegularOctahedron::orientedVertices;
std::array<Vector<3>, 4> RegularOctahedron::orientedFaceAxes;
std::array<Vector<3>, 6> RegularOctahedron::orientedEdgeAxes;

void RegularOctahedron::calculateStatic(const std::string &attr) {
    orientedVertices = {{Vector<3>{{ 1,  0,  0}} * edgeFactor,
                         Vector<3>{{-1,  0,  0}} * edgeFactor,
                         Vector<3>{{ 0,  1,  0}} * edgeFactor,
                         Vector<3>{{ 0, -1,  0}} * edgeFactor,
                         Vector<3>{{ 0,  0,  1}} * edgeFactor,
                         Vector<3>{{ 0,  0, -1}} * edgeFactor}};

    orientedEdgeAxes = {{orientedVertices[4] - orientedVertices[2],
                         orientedVertices[4] - orientedVertices[1],
                         orientedVertices[3] - orientedVertices[4],
                         orientedVertices[0] - orientedVertices[4],
                         orientedVertices[0] - orientedVertices[3],
                         orientedVertices[2] - orientedVertices[0]}};

    orientedFaceAxes = {{orientedEdgeAxes[3] ^ orientedEdgeAxes[0],
                         orientedEdgeAxes[0] ^ orientedEdgeAxes[1],
                         orientedEdgeAxes[2] ^ orientedEdgeAxes[1],
                         orientedEdgeAxes[2] ^ orientedEdgeAxes[3]}};
}

double RegularOctahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xHalfsize = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yHalfsize = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zHalfsize = std::abs(this->getOrientationMatrix().column(2) * axis);

    return std::max(std::max(xHalfsize, yHalfsize), zHalfsize) * edgeFactor;
}

std::array<Vector<3>, 4> RegularOctahedron::getFaceAxes() const {
    return this->faceAxes;
}

std::array<Vector<3>, 6> RegularOctahedron::getEdgeAxes() const {
    return this->edgeAxes;
}

std::array<Vector<3>, 8> RegularOctahedron::getVertices() const {
    return this->vertices;
}

RegularOctahedron::RegularOctahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularOctahedron>{orientation} {
    this->calculateVerticesAndAxes();
}

intersection::polyhedron RegularOctahedron::getTriangles() const {
    return intersection::polyhedron{
            {{this->vertices[0], this->vertices[2], this->vertices[4]}},
            {{this->vertices[1], this->vertices[5], this->vertices[3]}},
            {{this->vertices[4], this->vertices[2], this->vertices[1]}},
            {{this->vertices[3], this->vertices[5], this->vertices[0]}},
            {{this->vertices[4], this->vertices[1], this->vertices[3]}},
            {{this->vertices[2], this->vertices[0], this->vertices[5]}},
            {{this->vertices[3], this->vertices[0], this->vertices[4]}},
            {{this->vertices[1], this->vertices[2], this->vertices[5]}}
    };
}
