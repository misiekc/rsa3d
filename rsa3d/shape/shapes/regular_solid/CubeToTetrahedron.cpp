//
// Created by PKua on 26.09.18.
//

#include "CubeToTetrahedron.h"
#include "Tetrahedron.h"
#include "Octahedron.h"
#include "../../../utils/Assertions.h"

namespace {
    const double EPSILON = 1e-10;

    inline bool approx_equal(double a, double b) {
        return std::abs(a - b) < EPSILON;
    }
}

const UnoptimizedSATOverlapRS CubeToTetrahedron::overlapStrategy{};

RegularSolidBase::ShapeData CubeToTetrahedron::calculateStatic(const std::string &attr) {
    double aa, ac;
    std::istringstream attrStream(attr);
    attrStream >> aa >> ac;
    ValidateMsg(attrStream, "Expected 2 cut parameters");
    ValidateMsg(aa >= 0 && aa <= 1, "aa should be in [0, 1] range");
    ValidateMsg(ac >= 0 && ac <= 1, "ac should be in [0, 1] range");
    
    if (ac < aa)    std::swap(ac, aa);
    
    if (approx_equal(aa, 1) && approx_equal(ac, 1))
        throw ValidationException("aa = 1, ac = 1; this class doesn't support cubes; use Cube instead");
    else if (approx_equal(ac, 1) && approx_equal(aa, 0))
        return Tetrahedron::calculateStatic("dual");   // Use a dual to preserve continuity and be backwards compatible
    else if (approx_equal(aa, 0) && approx_equal(ac, 0))
        return Octahedron::calculateStatic("");
    else if (approx_equal(aa, 0))
        return calculateTruncatedTetrahedron(ac);
    else if (approx_equal(ac, 1))
        return calculateTruncatedCube(aa);
    else if (approx_equal(aa + ac, 1))
        return calculateCuboctahedron(aa);
    else if (aa + ac > 1)
        return calculateNonintersectingTruncations(aa, ac);
    else
        return calculateIntersectingTruncations(aa, ac);
}

OverlapStrategy<3, 0> *CubeToTetrahedron::createStrategy(const std::string &name) const {
    if (name == "sat")
        return new UnoptimizedSATOverlapRS();
    else
        return RegularSolid<CubeToTetrahedron>::createStrategy(name);
}

RegularSolidBase::ShapeData CubeToTetrahedron::calculateTruncatedTetrahedron(double ac) {
    ShapeData shapeData;

    shapeData.orientedVertices =
            {{{-1, -ac, -ac}}, {{-ac, -1, -ac}}, {{-ac, -ac, -1}}, {{-1,  ac, ac}}, {{ac, -1,  ac}}, {{ ac, ac, -1}},
             {{ 1,  ac, -ac}}, {{-ac,  1,  ac}}, {{ ac, -ac,  1}}, {{ 1, -ac, ac}}, {{ac,  1, -ac}}, {{-ac, ac,  1}}};

    shapeData.orientedFaces =
            {{10, 7, 3, 0, 2, 5}, {1, 4, 9, 6, 5, 2}, {10, 6, 9, 8, 11, 7}, {8, 4, 1, 0, 3, 11},

             {0, 1, 2}, {5, 6, 10}, {11, 3, 7}, {8, 9, 4}};

    return shapeData;
}

RegularSolidBase::ShapeData CubeToTetrahedron::calculateTruncatedCube(double aa) {
    double ba = 2*aa - 1;   // 0 means full truncation, 1 means none

    ShapeData shapeData;

    shapeData.orientedVertices =
            {{{ba,  1,  1}}, {{-ba,  -1,  1}}, {{-ba,  1,  -1}}, {{ba,  -1,  -1}},
             {{ 1, ba,  1}}, {{ -1, -ba,  1}}, {{ -1, ba,  -1}}, {{ 1, -ba,  -1}},
             {{ 1,  1, ba}}, {{ -1,  -1, ba}}, {{ -1,  1, -ba}}, {{ 1,  -1, -ba}},

             {{1, 1, -1}}, {{1, -1, 1}}, {{-1, 1, 1}}, {{-1, -1, -1}}};

    shapeData.orientedFaces =
            {{1, 13, 4, 0, 14, 5}, { 8, 12,  2, 10, 14, 0}, {14, 10, 6, 15,  9,  5},
             {2, 12, 7, 3, 15, 6}, {15,  3, 11, 13,  1, 9}, { 7, 12, 8,  4, 13, 11},

             {1, 5, 9}, {4, 8, 0}, {10, 2, 6}, {3, 7, 11}};

    return shapeData;
}

RegularSolidBase::ShapeData CubeToTetrahedron::calculateCuboctahedron(double aa) {
    double ba = 2*aa - 1;   // 0 means full truncation, 1 means none

    ShapeData shapeData;

    shapeData.orientedVertices =
            {{{ ba,  1, 1}}, {{1,  ba,  1}}, {{ 1, 1,  ba}}, {{ ba, -1, -1}}, {{-1,  ba, -1}}, {{-1, -1,  ba}},
             {{-ba, -1, 1}}, {{1, -ba, -1}}, {{-1, 1, -ba}}, {{-ba,  1, -1}}, {{-1, -ba,  1}}, {{ 1, -1, -ba}}};

    shapeData.orientedFaces =
            {{10, 8, 4, 5}, {6, 1, 0, 10}, {0, 2, 9, 8}, {3, 11, 6, 5}, {11, 7, 2, 1}, {7, 3, 4, 9},
             {6, 10, 5}, {10, 0, 8}, {8, 9, 4}, {4, 3, 5}, {11, 1, 6}, {3, 7, 11}, {1, 2, 0}, {7, 9, 2}};

    return shapeData;
}

RegularSolidBase::ShapeData CubeToTetrahedron::calculateNonintersectingTruncations(double aa, double ac) {
    double bc = 2*ac - 1;   // 0 means full truncation, 1 means none
    double ba = 2*aa - 1;

    ShapeData shapeData;

    shapeData.orientedVertices =
            {{{ba,  1,  1}}, {{-ba,  -1,  1}}, {{-ba,  1,  -1}}, {{ba,  -1,  -1}},
             {{ 1, ba,  1}}, {{ -1, -ba,  1}}, {{ -1, ba,  -1}}, {{ 1, -ba,  -1}},
             {{ 1,  1, ba}}, {{ -1,  -1, ba}}, {{ -1,  1, -ba}}, {{ 1,  -1, -ba}},
                    
             {{-bc,  -1,  -1}}, {{bc,  1,  -1}}, {{bc,  -1,  1}}, {{-bc,  1,  1}},
             {{ -1, -bc,  -1}}, {{ 1, bc,  -1}}, {{ 1, -bc,  1}}, {{ -1, bc,  1}},
             {{ -1,  -1, -bc}}, {{ 1,  1, -bc}}, {{ 1,  -1, bc}}, {{ -1,  1, bc}}};

    shapeData.orientedFaces =
            {{14, 18, 4,  0, 15, 19,  5,  1}, { 0,  8, 21, 13,  2, 10, 23, 15},
             {22, 11, 7, 17, 21,  8,  4, 18}, { 3, 12, 16,  6,  2, 13, 17,  7},
             {16, 20, 9,  5, 19, 23, 10,  6}, {12,  3, 11, 22, 14,  1,  9, 20},

             {16, 12, 20}, {9, 1, 5}, {3, 7, 11}, {22, 18, 14}, {17, 13, 21}, {2, 6, 10}, {23, 19, 15}, {8, 0, 4}};

    return shapeData;
}

RegularSolidBase::ShapeData CubeToTetrahedron::calculateIntersectingTruncations(double aa, double ac) {
    double p = aa + ac;
    double m = ac - aa;

    ShapeData shapeData;

    shapeData.orientedVertices =
            {{{-1, -p, -m}}, {{-1, -m, -p}}, {{-p, -m, -1}}, {{-m, -p, -1}}, {{-m, -1, -p}}, {{-p, -1, -m}},
             {{-1,  p,  m}}, {{-1,  m,  p}}, {{ p,  m, -1}}, {{ m,  p, -1}}, {{ m, -1,  p}}, {{ p, -1,  m}},
             {{ 1, -p,  m}}, {{ 1,  m, -p}}, {{-p,  m,  1}}, {{ m, -p,  1}}, {{ m,  1, -p}}, {{-p,  1,  m}},
             {{ 1,  p, -m}}, {{ 1, -m,  p}}, {{ p, -m,  1}}, {{-m,  p,  1}}, {{-m,  1,  p}}, {{ p,  1, -m}}};

    shapeData.orientedFaces =
            {{ 5,  4,  3,  2, 1, 0}, { 0,  7, 14, 15, 10, 5}, {11, 12, 13, 8, 3,  4}, {15, 20, 19, 12, 11, 10},
             {13, 18, 23, 16, 9, 8}, {16, 17,  6,  1,  2, 9}, {22, 21, 14, 7, 6, 17}, {23, 18, 19, 20, 21, 22},

             {19, 18, 13, 12}, {5, 10, 11, 4}, {8, 9, 2, 3}, {6, 7, 0, 1}, {23, 22, 17, 16}, {21, 20, 15, 14}};

    return shapeData;
}
