//
// Created by PKua on 21.04.18.
//

#include "RegularIcosahedron.h"

void RegularIcosahedron::calculateStatic(const std::string &attr) {
    auto &vertices = PlatonicSolid<RegularIcosahedron>::orientedVertices;
    auto &edgeAxes = PlatonicSolid<RegularIcosahedron>::orientedEdgeAxes;
    auto &faceAxes = PlatonicSolid<RegularIcosahedron>::orientedFaceAxes;
    
    vertices = {Vector<3>{{0,  1,  gold}} * edgeFactor,
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
                         Vector<3>{{ 1, -gold, 0}} * edgeFactor};

    edgeAxes = {vertices[0] - vertices[8],
                         vertices[0] - vertices[9],
                         vertices[0] - vertices[7],
                         vertices[0] - vertices[1],
                         vertices[0] - vertices[4],
                         vertices[4] - vertices[1],
                         vertices[8] - vertices[4],
                         vertices[9] - vertices[8],
                         vertices[7] - vertices[9],
                         vertices[1] - vertices[7],
                         vertices[1] - vertices[11],
                         vertices[4] - vertices[5],
                         vertices[8] - vertices[3],
                         vertices[9] - vertices[6],
                         vertices[7] - vertices[10]};

    faceAxes = {edgeAxes[5] ^ edgeAxes[3],
                         edgeAxes[6] ^ edgeAxes[4],
                         edgeAxes[7] ^ edgeAxes[1],
                         edgeAxes[8] ^ edgeAxes[2],
                         edgeAxes[9] ^ edgeAxes[3],
                         edgeAxes[5] ^ edgeAxes[10],
                         edgeAxes[6] ^ edgeAxes[11],
                         edgeAxes[7] ^ edgeAxes[12],
                         edgeAxes[8] ^ edgeAxes[13],
                         edgeAxes[9] ^ edgeAxes[14]};
}

RegularIcosahedron::RegularIcosahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularIcosahedron>{orientation} {
    this->calculateVerticesAndAxes();
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

intersection::tri_polyh RegularIcosahedron::getTriangles() const {
    auto vertices = this->getVertices();
    return intersection::tri_polyh{
            {{vertices[4], vertices[0], vertices[1]}},
            {{vertices[8], vertices[0], vertices[4]}},
            {{vertices[9], vertices[0], vertices[8]}},
            {{vertices[7], vertices[0], vertices[9]}},
            {{vertices[1], vertices[0], vertices[7]}},
            {{vertices[2], vertices[6], vertices[3]}},
            {{vertices[2], vertices[10], vertices[6]}},
            {{vertices[2], vertices[11], vertices[10]}},
            {{vertices[2], vertices[5], vertices[11]}},
            {{vertices[2], vertices[3], vertices[5]}},
            {{vertices[11], vertices[4], vertices[1]}},
            {{vertices[11], vertices[5], vertices[4]}},
            {{vertices[5], vertices[8], vertices[4]}},
            {{vertices[5], vertices[3], vertices[8]}},
            {{vertices[3], vertices[9], vertices[8]}},
            {{vertices[3], vertices[6], vertices[9]}},
            {{vertices[6], vertices[7], vertices[9]}},
            {{vertices[6], vertices[10], vertices[7]}},
            {{vertices[10], vertices[1], vertices[7]}},
            {{vertices[10], vertices[11], vertices[1]}}
    };
}