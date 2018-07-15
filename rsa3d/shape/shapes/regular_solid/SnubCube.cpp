//
// Created by PKua on 12.07.18.
//

#include "SnubCube.h"
#include "UnoptimizedSATOverlap.h"

void SnubCube::calculateStatic(const std::string &attr) {
    const double t = SnubCube::tribonacciConstant;
    const double u = 1/t;

    RegularSolid<SnubCube>::orientedVertices =
            {{{1, u,  t}}, {{1, -u, -t}}, {{-1, u, -t}}, {{-1, -u,  t}},
             {{1, t, -u}}, {{1, -t,  u}}, {{-1, t,  u}}, {{-1, -t, -u}},

             {{u, t,  1}}, {{u, -t, -1}}, {{-u, t, -1}}, {{-u, -t,  1}},
             {{t, u, -1}}, {{t, -u,  1}}, {{-t, u,  1}}, {{-t, -u, -1}},

             {{t, 1,  u}}, {{t, -1, -u}}, {{-t, 1, -u}}, {{-t, -1,  u}},
             {{u, 1, -t}}, {{u, -1,  t}}, {{-u, 1,  t}}, {{-u, -1, -t}}};

    RegularSolid<SnubCube>::orientedFaces =
            {{0, 22, 3, 21}, {8, 4, 10, 6}, {20, 1, 23, 2}, {9, 5, 11, 7}, {15, 19, 14, 18}, {13, 17, 12, 16},

             {12, 20,  4}, {12,  4, 16}, {16,  4,  8}, {16,  8,  0}, {13, 16,  0}, {13,  0, 21}, { 5, 13, 21}, { 5, 17, 13},
             { 9, 17,  5}, { 9,  1, 17}, { 1, 12, 17}, { 1, 20, 12}, {20, 10,  4}, {20,  2, 10}, { 6, 22,  8}, { 8, 22,  0},
             { 3, 11, 21}, {21, 11,  5}, { 7, 23,  9}, { 9, 23,  1}, { 7, 19, 15}, {11, 19,  7}, {11,  3, 19}, { 3, 14, 19},
             { 3, 22, 14}, {22,  6, 14}, {14,  6, 18}, { 6, 10, 18}, {18, 10,  2}, {18,  2, 15}, {15,  2, 23}, {15, 23,  7}};
}

double SnubCube::projectionHalfsize(const Vector<3> &axis) const {
    throw std::runtime_error("unimplemented");
}

OverlapStrategy<3, 0> *SnubCube::createStrategy(const std::string &name) const {
    if (name == "sat")
        return new UnoptimizedSATOverlap<SnubCube>();
    else
        return RegularSolid<SnubCube>::createStrategy(name);
}

bool SnubCube::overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const {
    SnubCube other = dynamic_cast<const SnubCube&>(*s);     // Make a copy
    this->applyBC(bc, &other);
    UnoptimizedSATOverlap<SnubCube> satOverlap;
    return satOverlap.overlap(this, &other);
}
