//
// Created by PKua on 03.12.17.
//

#include "TriTriOverlapCB.h"


bool TriTriOverlapCB::overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const {
    auto cube1 = dynamic_cast<const Cuboid*>(first);
    auto cube2 = dynamic_cast<const Cuboid*>(second);

    return collision::polyh_polyh(cube1->obtainTris(), cube2->obtainTris());
}

std::string TriTriOverlapCB::getName() const {
    return "TriTriOverlapCB";
}

void TriTriOverlapCB::runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) const {
    cube1->obtainTris();
    cube2->obtainTris();
}
