//
// Created by pkua on 17.03.19.
//

#include "Stolen2DOverlapSC.h"
#include "Spherocylinder.h"
#include "../../../boundary_conditions/FreeBC.h"

SpheroCylinder2D toSpheroCylinder2D(const Spherocylinder<2> &sc) {
    Orientation<1> rotation = {{sc.getAngle()}};
    Vector<2> position = sc.getPosition();

    SpheroCylinder2D result;
    result.translate(position);
    result.rotate(rotation);
    return result;
}

bool Stolen2DOverlapSC::overlap(const Shape<2, 0> *first, const Shape<2, 0> *second) const {
    auto firstSC = toSpheroCylinder2D(dynamic_cast<const Spherocylinder<2> &>(*first));
    auto secondSC = toSpheroCylinder2D(dynamic_cast<const Spherocylinder<2> &>(*second));

    FreeBC<2> bc;
    return firstSC.overlap(&bc, &secondSC);
}
