//
// Created by PKua on 03.06.18.
//

#ifndef RSA3D_ISOTROPICFACTORY_H
#define RSA3D_ISOTROPICFACTORY_H


#include "ShapePairFactory.h"
#include "Distrubution.h"


class ParallelPairFactory : public RSAShapePairFactory {
private:
    Distribution &distribution;

public:
    explicit ParallelPairFactory(Distribution &distribution);

    ShapePair generate() override;
    std::string getDescription() const override;
};

#endif //RSA3D_ISOTROPICFACTORY_H
