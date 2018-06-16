//--------------------------------------------------------------------------------------------
// Factory generating cuboids from ball of given radius
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _BALL_FACTORY_H
    #define _BALL_FACTORY_H

#include "ShapePairFactory.h"

/**
 * @brief ShapePairFactory generating shapes from a hyperball of given radius.
 *
 * A dimension is determined by current @a RSA_SPATIAL_DIMESION.
 */
class BallFactory : public RSAShapePairFactory
{
private:
    double radius{0.5};
    unsigned int no{};

    /* Helper method. Creates random Cuboid based on objects parameters. Delegate shape creation to the standard
     * ShapeFactory */
    std::unique_ptr<RSAShape> randomShape();
    double randomGaussian();

public:
    ShapePair generate() override;
    std::string getDescription() const override;
    void setRadius(double _radius);
};

#endif // _BALL_FACTORY_H
