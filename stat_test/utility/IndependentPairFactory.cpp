//--------------------------------------------------------------------------------------------
// Factory generating cuboids from ball of given radius
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include <sstream>

#include "IndependentPairFactory.h"
#include "ShapeGenerators.h"

RSAShapePairFactory::ShapePair IndependentPairFactory::generate() {
    return {generate_randomly_oriented_shape(this->distribution.randomVector(&rnd), &rnd),
            generate_randomly_oriented_shape(this->distribution.randomVector(&rnd), &rnd)};
}

std::string IndependentPairFactory::getDescription() const {
    return "IndependentPairFactory over " + this->distribution.getDescription();
}

IndependentPairFactory::IndependentPairFactory(Distribution &distribution) : distribution(distribution) {}

Distribution &IndependentPairFactory::getDistribution() {
    return this->distribution;
}
