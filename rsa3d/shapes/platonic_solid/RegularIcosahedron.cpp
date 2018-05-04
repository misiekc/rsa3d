//
// Created by PKua on 21.04.18.
//

#include "RegularIcosahedron.h"

std::array<Vector<3>, 20> RegularIcosahedron::orientedVertices;
std::array<Vector<3>, 10> RegularIcosahedron::orientedFaceAxes;
std::array<Vector<3>, 15> RegularIcosahedron::orientedEdgeAxes;

void RegularIcosahedron::calculateStatic(const std::string &attr) {
    orientedVertices = {{Vector<3>{{0,  1,  gold}} * edgeFactor,
                         Vector<3>{{0, -1,  gold}} * edgeFactor,
                         Vector<3>{{0, -1, -gold}} * edgeFactor,
                         Vector<3>{{0,  1, -gold}} * edgeFactor,

                         Vector<3>{{ gold, 0,  1}} * edgeFactor,
                         Vector<3>{{ gold, 0, -1}} * edgeFactor,
                         Vector<3>{{-gold, 0, -1}} * edgeFactor,
                         Vector<3>{{-gold, 0,  1}} * edgeFactor,

                         Vector<3>{{ 1,  gold, 0}} * edgeFactor,
                         Vector<3>{{-1,  gold, 0}} * edgeFactor,
                         Vector<3>{{-1, -gold, 0}} * edgeFactor,
                         Vector<3>{{ 1, -gold, 0}} * edgeFactor}};

    orientedEdgeAxes = {{orientedVertices[0] - orientedVertices[8],
                         orientedVertices[0] - orientedVertices[9],
                         orientedVertices[0] - orientedVertices[7],
                         orientedVertices[0] - orientedVertices[1],
                         orientedVertices[0] - orientedVertices[4],
                         orientedVertices[4] - orientedVertices[1],
                         orientedVertices[8] - orientedVertices[4],
                         orientedVertices[9] - orientedVertices[8],
                         orientedVertices[7] - orientedVertices[9],
                         orientedVertices[1] - orientedVertices[7],
                         orientedVertices[1] - orientedVertices[11],
                         orientedVertices[4] - orientedVertices[5],
                         orientedVertices[8] - orientedVertices[3],
                         orientedVertices[9] - orientedVertices[6],
                         orientedVertices[7] - orientedVertices[10]}};

    orientedFaceAxes = {{orientedEdgeAxes[5] ^ orientedEdgeAxes[3],
                         orientedEdgeAxes[6] ^ orientedEdgeAxes[4],
                         orientedEdgeAxes[7] ^ orientedEdgeAxes[1],
                         orientedEdgeAxes[8] ^ orientedEdgeAxes[2],
                         orientedEdgeAxes[9] ^ orientedEdgeAxes[3],
                         orientedEdgeAxes[5] ^ orientedEdgeAxes[10],
                         orientedEdgeAxes[6] ^ orientedEdgeAxes[11],
                         orientedEdgeAxes[7] ^ orientedEdgeAxes[12],
                         orientedEdgeAxes[8] ^ orientedEdgeAxes[13],
                         orientedEdgeAxes[9] ^ orientedEdgeAxes[14]}};
}

RegularIcosahedron::RegularIcosahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularIcosahedron>{orientation} {
    this->calculateVerticesAndAxes();
}

std::array<Vector<3>, 20> RegularIcosahedron::getVertices() const {
    return this->vertices;
}

std::array<Vector<3>, 10> RegularIcosahedron::getFaceAxes() const {
    return this->faceAxes;
}

std::array<Vector<3>, 15> RegularIcosahedron::getEdgeAxes() const {
    return this->edgeAxes;
}

double RegularIcosahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xHalfsize = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yHalfsize = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zHalfsize = std::abs(this->getOrientationMatrix().column(2) * axis);

    double xRectHalfsize = yHalfsize + zHalfsize * gold;
    double yRectHalfsize = zHalfsize + xHalfsize * gold;
    double zRectHalfsize = xHalfsize + yHalfsize * gold;

    return std::max(std::max(xRectHalfsize, yRectHalfsize), zRectHalfsize) * edgeFactor;
}

intersection::polyhedron RegularIcosahedron::getTriangles() const {
    return intersection::polyhedron{
            {{this->vertices[4], this->vertices[0], this->vertices[1]}},
            {{this->vertices[8], this->vertices[0], this->vertices[4]}},
            {{this->vertices[9], this->vertices[0], this->vertices[8]}},
            {{this->vertices[7], this->vertices[0], this->vertices[9]}},
            {{this->vertices[1], this->vertices[0], this->vertices[7]}},
            {{this->vertices[2], this->vertices[6], this->vertices[3]}},
            {{this->vertices[2], this->vertices[10], this->vertices[6]}},
            {{this->vertices[2], this->vertices[11], this->vertices[10]}},
            {{this->vertices[2], this->vertices[5], this->vertices[11]}},
            {{this->vertices[2], this->vertices[3], this->vertices[5]}},
            {{this->vertices[11], this->vertices[4], this->vertices[1]}},
            {{this->vertices[11], this->vertices[5], this->vertices[4]}},
            {{this->vertices[5], this->vertices[8], this->vertices[4]}},
            {{this->vertices[5], this->vertices[3], this->vertices[8]}},
            {{this->vertices[3], this->vertices[9], this->vertices[8]}},
            {{this->vertices[3], this->vertices[6], this->vertices[9]}},
            {{this->vertices[6], this->vertices[7], this->vertices[9]}},
            {{this->vertices[6], this->vertices[10], this->vertices[7]}},
            {{this->vertices[10], this->vertices[1], this->vertices[7]}},
            {{this->vertices[10], this->vertices[11], this->vertices[1]}}
    };
}