//
// Created by PKua on 11.08.18.
//

#include "ParallelPairFactory.h"
#include "ShapeGenerators.h"

ParallelPairFactory::ParallelPairFactory(Distribution &distribution) : distribution(distribution) {}

RSAShapePairFactory::ShapePair ParallelPairFactory::generate() {
    auto randomShape = generate_randomly_oriented_shape(this->distribution.randomVector(&rnd), &rnd);

    auto shape1 = randomShape->clone();
    auto shape2 = randomShape->clone();
    shape2->translate(this->distribution.randomVector(&rnd) - randomShape->getPosition());

    return ShapePair(shape1, shape2);
}

std::string ParallelPairFactory::getDescription() const {
    return "ParallelPairFactory over " + distribution.getDescription();
}
