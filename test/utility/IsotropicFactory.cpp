//
// Created by PKua on 03.06.18.
//

#include "IsotropicFactory.h"
#include "../../rsa3d/Vector.h"

ShapePairFactory::ShapePair IsotropicFactory::generate() {
    auto pair = factory.generate();

    using V = Vector<RSA_SPATIAL_DIMENSION>;
    double transArr[RSA_SPATIAL_DIMENSION];
    V translation = V(pair.second()->getPosition()) - V(pair.first()->getPosition());
    translation.copyToArray(transArr);

    auto shape1 = pair.first()->clone();
    auto shape2 = pair.first()->clone();
    shape2->translate(transArr);
    return ShapePair(shape1, shape2);
}

std::string IsotropicFactory::getDescription() const {
    return "IsotropicFactory over " + factory.getDescription();
}
