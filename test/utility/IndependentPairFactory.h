//--------------------------------------------------------------------------------------------
// Factory generating cuboids from ball of given radius
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _BALL_FACTORY_H
    #define _BALL_FACTORY_H

#include "ShapePairFactory.h"
#include "Distrubution.h"


class IndependentPairFactory : public RSAShapePairFactory
{
private:
    Distribution &distribution;

public:
    explicit IndependentPairFactory(Distribution &distribution);

    ShapePair generate() override;
    std::string getDescription() const override;

    Distribution &getDistribution();
};

#endif // _BALL_FACTORY_H
