//
// Created by PKua on 21.04.18.
//

#include "RegularDodecahedron.h"

std::array<Vector<3>, 20> RegularDodecahedron::orientedVertices;
std::array<Vector<3>, 6> RegularDodecahedron::orientedFaceAxes;
std::array<Vector<3>, 15> RegularDodecahedron::orientedEdgeAxes;

void RegularDodecahedron::calculateStatic(const std::string &attr) {
    orientedVertices = {{Vector<3>{{ 1,  1,  1}} * edgeFactor,
                         Vector<3>{{-1,  1,  1}} * edgeFactor,
                         Vector<3>{{-1, -1,  1}} * edgeFactor,
                         Vector<3>{{ 1, -1,  1}} * edgeFactor,
                         Vector<3>{{-1, -1, -1}} * edgeFactor,
                         Vector<3>{{-1,  1, -1}} * edgeFactor,
                         Vector<3>{{ 1,  1, -1}} * edgeFactor,
                         Vector<3>{{ 1, -1, -1}} * edgeFactor,

                         Vector<3>{{ gold,  1/gold, 0}} * edgeFactor,
                         Vector<3>{{-gold,  1/gold, 0}} * edgeFactor,
                         Vector<3>{{-gold, -1/gold, 0}} * edgeFactor,
                         Vector<3>{{ gold, -1/gold, 0}} * edgeFactor,

                         Vector<3>{{0,  gold,  1/gold}} * edgeFactor,
                         Vector<3>{{0, -gold,  1/gold}} * edgeFactor,
                         Vector<3>{{0, -gold, -1/gold}} * edgeFactor,
                         Vector<3>{{0,  gold, -1/gold}} * edgeFactor,

                         Vector<3>{{ 1/gold, 0,  gold}} * edgeFactor,
                         Vector<3>{{ 1/gold, 0, -gold}} * edgeFactor,
                         Vector<3>{{-1/gold, 0, -gold}} * edgeFactor,
                         Vector<3>{{-1/gold, 0,  gold}} * edgeFactor}};

    orientedEdgeAxes = {{orientedVertices[0] - orientedVertices[8],
                         orientedVertices[0] - orientedVertices[16],
                         orientedVertices[16] - orientedVertices[3],
                         orientedVertices[3] - orientedVertices[11],
                         orientedVertices[8] - orientedVertices[11],
                         orientedVertices[0] - orientedVertices[12],
                         orientedVertices[16] - orientedVertices[19],
                         orientedVertices[3] - orientedVertices[13],
                         orientedVertices[11] - orientedVertices[7],
                         orientedVertices[8] - orientedVertices[6],
                         orientedVertices[7] - orientedVertices[14],
                         orientedVertices[2] - orientedVertices[13],
                         orientedVertices[19] - orientedVertices[1],
                         orientedVertices[12] - orientedVertices[15],
                         orientedVertices[6] - orientedVertices[17]}};

    orientedFaceAxes = {{orientedEdgeAxes[1] ^ orientedEdgeAxes[0],
                         orientedEdgeAxes[14] ^ orientedEdgeAxes[9],
                         orientedEdgeAxes[0] ^ orientedEdgeAxes[5],
                         orientedEdgeAxes[5] ^ orientedEdgeAxes[1],
                         orientedEdgeAxes[6] ^ orientedEdgeAxes[2],
                         orientedEdgeAxes[7] ^ orientedEdgeAxes[3]}};
}

RegularDodecahedron::RegularDodecahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularDodecahedron>{orientation} {
    this->calculateVerticesAndAxes();
}

std::array<Vector<3>, 20> RegularDodecahedron::getVertices() const {
    return this->vertices;
}

std::array<Vector<3>, 6> RegularDodecahedron::getFaceAxes() const {
    return this->faceAxes;
}

std::array<Vector<3>, 15> RegularDodecahedron::getEdgeAxes() const {
    return this->edgeAxes;
}

double RegularDodecahedron::projectionHalfsize(const Vector<3> &axis) const {
    double xHalfsize = std::abs(this->getOrientationMatrix().column(0) * axis);
    double yHalfsize = std::abs(this->getOrientationMatrix().column(1) * axis);
    double zHalfsize = std::abs(this->getOrientationMatrix().column(2) * axis);

    double cubeHalfsize = xHalfsize + yHalfsize + zHalfsize;
    double xRectHalfsize = xHalfsize * gold + yHalfsize / gold;
    double yRectHalfsize = yHalfsize * gold + zHalfsize / gold;
    double zRectHalfsize = zHalfsize * gold + xHalfsize / gold;

    return std::max(std::max(std::max(cubeHalfsize, xRectHalfsize), yRectHalfsize), zRectHalfsize) * edgeFactor;
}

intersection::polyhedron RegularDodecahedron::getTriangles() const {
    return intersection::polyhedron{
            {{this->vertices[3], this->vertices[0], this->vertices[16]}},
            {{this->vertices[3], this->vertices[11], this->vertices[0]}},
            {{this->vertices[11], this->vertices[8], this->vertices[0]}},

            {{this->vertices[17], this->vertices[6], this->vertices[7]}},
            {{this->vertices[7], this->vertices[6], this->vertices[8]}},
            {{this->vertices[7], this->vertices[8], this->vertices[11]}},

            {{this->vertices[15], this->vertices[12], this->vertices[0]}},
            {{this->vertices[15], this->vertices[0], this->vertices[8]}},
            {{this->vertices[15], this->vertices[8], this->vertices[6]}},

            {{this->vertices[16], this->vertices[0], this->vertices[12]}},
            {{this->vertices[16], this->vertices[12], this->vertices[1]}},
            {{this->vertices[16], this->vertices[1], this->vertices[19]}},

            {{this->vertices[3], this->vertices[16], this->vertices[19]}},
            {{this->vertices[3], this->vertices[19], this->vertices[2]}},
            {{this->vertices[3], this->vertices[2], this->vertices[13]}},

            {{this->vertices[11], this->vertices[3], this->vertices[13]}},
            {{this->vertices[11], this->vertices[13], this->vertices[14]}},
            {{this->vertices[11], this->vertices[14], this->vertices[7]}},

            {{this->vertices[4], this->vertices[10], this->vertices[9]}},
            {{this->vertices[4], this->vertices[9], this->vertices[5]}},
            {{this->vertices[4], this->vertices[5], this->vertices[18]}},

            {{this->vertices[2], this->vertices[19], this->vertices[1]}},
            {{this->vertices[2], this->vertices[1], this->vertices[9]}},
            {{this->vertices[2], this->vertices[9], this->vertices[10]}},

            {{this->vertices[4], this->vertices[14], this->vertices[13]}},
            {{this->vertices[4], this->vertices[13], this->vertices[2]}},
            {{this->vertices[4], this->vertices[2], this->vertices[10]}},

            {{this->vertices[18], this->vertices[17], this->vertices[7]}},
            {{this->vertices[18], this->vertices[7], this->vertices[14]}},
            {{this->vertices[18], this->vertices[14], this->vertices[4]}},

            {{this->vertices[17], this->vertices[18], this->vertices[5]}},
            {{this->vertices[17], this->vertices[5], this->vertices[15]}},
            {{this->vertices[17], this->vertices[15], this->vertices[6]}},

            {{this->vertices[9], this->vertices[1], this->vertices[12]}},
            {{this->vertices[9], this->vertices[12], this->vertices[15]}},
            {{this->vertices[9], this->vertices[15], this->vertices[5]}}
    };
}

std::vector<double> RegularDodecahedron::calculateOrder(const OrderCalculable *other) const {
    auto &otherDod = dynamic_cast<const RegularDodecahedron&>(*other);

    double cos6sum = 0;
    double maxAbsCos = 0;
    for (const auto &a1 : this->faceAxes) {
        for (const auto &a2 : otherDod.faceAxes) {
            double absCos = std::abs(a1 * a2);
            cos6sum += std::pow(absCos, 6);
            if (absCos > maxAbsCos)
                maxAbsCos = absCos;
        }
    }
    double nematicOrder = P2(maxAbsCos);
    double dodecahedralOrder = 25./192 * (7*cos6sum - 36);

    return {nematicOrder, dodecahedralOrder};
}
