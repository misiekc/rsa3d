//
// Created by PKua on 21.04.18.
//

#include "RegularDodecahedron.h"

void RegularDodecahedron::calculateStatic(const std::string &attr) {
    auto &vertices = PlatonicSolid<RegularDodecahedron>::orientedVertices;
    auto &edgeAxes = PlatonicSolid<RegularDodecahedron>::orientedEdgeAxes;
    auto &faceAxes = PlatonicSolid<RegularDodecahedron>::orientedFaceAxes;
    
    vertices = {Vector<3>{{ 1,  1,  1}} * edgeFactor,
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
                         Vector<3>{{-1/gold, 0,  gold}} * edgeFactor};

    edgeAxes = {vertices[0] - vertices[8],
                         vertices[0] - vertices[16],
                         vertices[16] - vertices[3],
                         vertices[3] - vertices[11],
                         vertices[8] - vertices[11],
                         vertices[0] - vertices[12],
                         vertices[16] - vertices[19],
                         vertices[3] - vertices[13],
                         vertices[11] - vertices[7],
                         vertices[8] - vertices[6],
                         vertices[7] - vertices[14],
                         vertices[2] - vertices[13],
                         vertices[19] - vertices[1],
                         vertices[12] - vertices[15],
                         vertices[6] - vertices[17]};

    faceAxes = {edgeAxes[1] ^ edgeAxes[0],
                         edgeAxes[14] ^ edgeAxes[9],
                         edgeAxes[0] ^ edgeAxes[5],
                         edgeAxes[5] ^ edgeAxes[1],
                         edgeAxes[6] ^ edgeAxes[2],
                         edgeAxes[7] ^ edgeAxes[3]};
}

RegularDodecahedron::RegularDodecahedron(const Matrix<3, 3> &orientation) : PlatonicSolid<RegularDodecahedron>{orientation} {
    this->calculateVerticesAndAxes();
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

intersection::tri_polyh RegularDodecahedron::getTriangles() const {
    auto vertices = this->getVertices();
    return intersection::tri_polyh{
            {{vertices[3], vertices[0], vertices[16]}},
            {{vertices[3], vertices[11], vertices[0]}},
            {{vertices[11], vertices[8], vertices[0]}},

            {{vertices[17], vertices[6], vertices[7]}},
            {{vertices[7], vertices[6], vertices[8]}},
            {{vertices[7], vertices[8], vertices[11]}},

            {{vertices[15], vertices[12], vertices[0]}},
            {{vertices[15], vertices[0], vertices[8]}},
            {{vertices[15], vertices[8], vertices[6]}},

            {{vertices[16], vertices[0], vertices[12]}},
            {{vertices[16], vertices[12], vertices[1]}},
            {{vertices[16], vertices[1], vertices[19]}},

            {{vertices[3], vertices[16], vertices[19]}},
            {{vertices[3], vertices[19], vertices[2]}},
            {{vertices[3], vertices[2], vertices[13]}},

            {{vertices[11], vertices[3], vertices[13]}},
            {{vertices[11], vertices[13], vertices[14]}},
            {{vertices[11], vertices[14], vertices[7]}},

            {{vertices[4], vertices[10], vertices[9]}},
            {{vertices[4], vertices[9], vertices[5]}},
            {{vertices[4], vertices[5], vertices[18]}},

            {{vertices[2], vertices[19], vertices[1]}},
            {{vertices[2], vertices[1], vertices[9]}},
            {{vertices[2], vertices[9], vertices[10]}},

            {{vertices[4], vertices[14], vertices[13]}},
            {{vertices[4], vertices[13], vertices[2]}},
            {{vertices[4], vertices[2], vertices[10]}},

            {{vertices[18], vertices[17], vertices[7]}},
            {{vertices[18], vertices[7], vertices[14]}},
            {{vertices[18], vertices[14], vertices[4]}},

            {{vertices[17], vertices[18], vertices[5]}},
            {{vertices[17], vertices[5], vertices[15]}},
            {{vertices[17], vertices[15], vertices[6]}},

            {{vertices[9], vertices[1], vertices[12]}},
            {{vertices[9], vertices[12], vertices[15]}},
            {{vertices[9], vertices[15], vertices[5]}}
    };
}

std::vector<double> RegularDodecahedron::calculateOrder(const OrderCalculable *other) const {
    auto &otherDod = dynamic_cast<const RegularDodecahedron&>(*other);

    double cos6sum = 0;
    double maxAbsCos = 0;
    auto faceAxes = this->getFaceAxes();
    auto otherFaceAxes = otherDod.getFaceAxes();
    for (const auto &a1 : faceAxes) {
        for (const auto &a2 : otherFaceAxes) {
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
