//
// Created by PKua on 26.09.18.
//

#include "CubeToTetrahedron.h"

const UnoptimizedSATOverlapRS CubeToTetrahedron::overlapStrategy{};

void CubeToTetrahedron::calculateStatic(const std::string &attr) {
    std::istringstream attrStream(attr);
    double xi;
    attrStream >> xi;
    if (!attrStream || xi < -1 || xi > 1)
        throw std::runtime_error("attr must contain xi parameter from [-1, 1]");

    RegularSolid<CubeToTetrahedron>::orientedVertices =
            {{{xi,  1,  1}}, {{-xi,  -1,  1}}, {{-xi,  1,  -1}}, {{xi,  -1,  -1}},
             {{ 1, xi,  1}}, {{ -1, -xi,  1}}, {{ -1, xi,  -1}}, {{ 1, -xi,  -1}},
             {{ 1,  1, xi}}, {{- 1,  -1, xi}}, {{ -1,  1, -xi}}, {{1 ,  -1, -xi}},

             {{1, 1, -1}}, {{1, -1, 1}}, {{-1, 1, 1}}, {{-1, -1, -1}}};

    RegularSolid<CubeToTetrahedron>::orientedFaces =
            {{1, 13, 4, 0, 14, 5}, { 8, 12,  2, 10, 14, 0}, {14, 10, 6, 15,  9,  5},
             {2, 12, 7, 3, 15, 6}, {15,  3, 11, 13,  1, 9}, { 7, 12, 8,  4, 13, 11},

             {1, 5, 9}, {4, 8, 0}, {10, 2, 6}, {3, 7, 11}};
}

OverlapStrategy<3, 0> *CubeToTetrahedron::createStrategy(const std::string &name) const {
    if (name == "sat")
        return new UnoptimizedSATOverlapRS();
    else
        return RegularSolid<CubeToTetrahedron>::createStrategy(name);
}