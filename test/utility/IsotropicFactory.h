//
// Created by PKua on 03.06.18.
//

#ifndef RSA3D_ISOTROPICFACTORY_H
#define RSA3D_ISOTROPICFACTORY_H


#include "ShapePairFactory.h"

class IsotropicFactory : public ShapePairFactory {
private:
    ShapePairFactory &factory;

public:
    IsotropicFactory(ShapePairFactory &factory) : factory(factory) {}

    ShapePair generate() override;
    std::string getDescription() const override;
};


#endif //RSA3D_ISOTROPICFACTORY_H
