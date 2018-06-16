//
// Created by PKua on 03.12.17.
//

#include "TriTriOverlap.h"


bool TriTriOverlap::overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const {
    auto cube1 = dynamic_cast<const Cuboid*>(first);
    auto cube2 = dynamic_cast<const Cuboid*>(second);

    return intersection::polyh_polyh(cube1->obtainTris(), cube2->obtainTris());
}

std::string TriTriOverlap::getName() const {
    return "TriTriOverlap";
}

void TriTriOverlap::runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) const {
    cube1->obtainTris();
    cube2->obtainTris();
}
