//
// Created by PKua on 03.06.18.
//

#ifndef RSA3D_ISOTROPICFACTORY_H
#define RSA3D_ISOTROPICFACTORY_H


#include "ShapePairFactory.h"
#include "../../rsa3d/Vector.h"

/**
 * @brief ShapePairFactory decorator ensuring both shapes have identical orientations.
 *
 * @tparam SD spatial dimension of produces shapes
 * @tparam AD angular dimesnion of produces shapes
 */
template <unsigned short SD, unsigned short AD>
class IsotropicFactory : public ShapePairFactory<SD, AD> {
private:
    using shape_pair_factory = ShapePairFactory<SD, AD>;
    using shape_pair = typename ShapePairFactory<SD, AD>::ShapePair;

    shape_pair_factory &factory;

public:
    IsotropicFactory(shape_pair_factory &factory) : factory(factory) {}

    shape_pair generate() override {
        auto pair = factory.generate();
        auto shape1 = pair.first()->clone();
        auto shape2 = pair.first()->clone();
        shape2->translate(pair.second()->getPosition() - pair.first()->getPosition());
        return shape_pair(shape1, shape2);
    };

    std::string getDescription() const override {
        return "IsotropicFactory over " + factory.getDescription();
    }
};

using RSAIsotropicFactory = IsotropicFactory<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>;

#endif //RSA3D_ISOTROPICFACTORY_H
