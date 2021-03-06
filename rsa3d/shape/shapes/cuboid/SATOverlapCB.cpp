//
// Created by PKua on 03.12.17.
//

#include "SATOverlapCB.h"


bool SATOverlapCB::overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const {
    auto cube1 = dynamic_cast<const Cuboid*>(first);
    auto cube2 = dynamic_cast<const Cuboid*>(second);

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
    cube2->obtainVertices(vertices2, Vector<3>());

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
bool SATOverlapCB::checkSeparatingAxis(const Vector<3> & _axis, Vector<3> * _vert1, Vector<3> * _vert2) const
{
    interval this_int = this->getProjection(_axis, _vert1);
    interval second_int = this->getProjection(_axis, _vert2);

    return std::min(this_int.second, second_int.second) >= std::max(this_int.first, second_int.first);
}

// Projects polyhedron _polyh on axis _axis and returns interval given by
// the projection
//----------------------------------------------------------------------------
SATOverlapCB::interval SATOverlapCB::getProjection(const Vector<3> & _axis, Vector<3> * _vert) const
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

std::string SATOverlapCB::getName() const {
    return "SATOverlapCB";
}

void SATOverlapCB::runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) const {
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
