//
// Created by PKua on 03.12.17.
//

#include "SATOverlap.h"


bool SATOverlap::overlap(const Cuboid *cube1, const Cuboid *cube2, BoundaryConditions *bc) {
    double trans_arr[3];
    Vector<3> translation(bc->getTranslation(trans_arr, cube1->getPosition(), cube2->getPosition()));
    Matrix<3, 3> orientation1 = cube1->getOrientation();
    Matrix<3, 3> orientation2 = cube2->getOrientation();
    Vector<3> vertices1[VERTEX::NUM_OF];
    Vector<3> vertices2[VERTEX::NUM_OF];

    double size[3];
    Cuboid::getSize(size);

    // Calculate axes orthogonal to separating plane
    Vector<3> axes1[] = {
            orientation1 * Vector<3>{{1, 0, 0}},
            orientation1 * Vector<3>{{0, 1, 0}},
            orientation1 * Vector<3>{{0, 0, 1}}
    };
    Vector<3> axes2[] = {
            orientation2 * Vector<3>{{1, 0, 0}},
            orientation2 * Vector<3>{{0, 1, 0}},
            orientation2 * Vector<3>{{0, 0, 1}}
    };

    cube1->obtainVertices(vertices1, Vector<3>());
    cube2->obtainVertices(vertices2, translation);

    // Check all possible separating axes - edge lines ...
    for (int i = 0; i < 3; i++)
        if (!this->checkSeparatingAxis (axes1[i], vertices1, vertices2))
            return false;
    for (int i = 0; i < 3; i++)
        if (!this->checkSeparatingAxis (axes2[i], vertices1, vertices2))
            return false;
    // ... and their cross products
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (!this->checkSeparatingAxis (axes1[i] ^ axes2[j], vertices1, vertices2))
                return false;
    return true;
}

// Checks whether this and _second projections on axis _axis overlap. If so,
// returns true
//----------------------------------------------------------------------------
bool SATOverlap::checkSeparatingAxis(const Vector<3> & _axis, Vector<3> * _vert1, Vector<3> * _vert2) const
{
    interval this_int = this->getProjection(_axis, _vert1);
    interval second_int = this->getProjection(_axis, _vert2);

    return std::min(this_int.second, second_int.second) >= std::max(this_int.first, second_int.first);
}

// Projects polyhedron _polyh on axis _axis and returns interval given by
// the projection
//----------------------------------------------------------------------------
SATOverlap::interval SATOverlap::getProjection(const Vector<3> & _axis, Vector<3> * _vert) const
{
    interval proj_int = {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};

    // Find enpoints of polyhedron projection (multiplied by unknown but const for _axit factor)
    double proj;
    for (unsigned long i = 0; i < VERTEX::NUM_OF; i++)
    {
        proj = _vert[i] * _axis;
        if (proj < proj_int.first)
            proj_int.first = proj;
        if (proj > proj_int.second)
            proj_int.second = proj;
    }
    return proj_int;
}

std::string SATOverlap::getName() {
    return "SATOverlap";
}

void SATOverlap::runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) {
    Matrix<3, 3> orientation1 = cube1->getOrientation();
    Matrix<3, 3> orientation2 = cube2->getOrientation();

    orientation1 * Vector<3>{{1, 0, 0}};
    orientation1 * Vector<3>{{0, 1, 0}};
    orientation1 * Vector<3>{{0, 0, 1}};

    orientation2 * Vector<3>{{1, 0, 0}};
    orientation2 * Vector<3>{{0, 1, 0}};
    orientation2 * Vector<3>{{0, 0, 1}};

    Vector<3> vertices1[VERTEX::NUM_OF];
    Vector<3> vertices2[VERTEX::NUM_OF];
    cube1->obtainVertices(vertices1, Vector<3>());
    cube2->obtainVertices(vertices2, Vector<3>());
}
