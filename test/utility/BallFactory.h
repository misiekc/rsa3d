//--------------------------------------------------------------------------------------------
// Factory generating cuboids from ball of given radius
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _BALL_FACTORY_H
    #define _BALL_FACTORY_H

#include "ShapePairFactory.h"

/**
 * @brief ShapePairFactory generating shapes from a ball of given size.
 */
class BallFactory : public ShapePairFactory
{
private:
    double radius{0.5};
    unsigned int no{};
    RND rnd{};

    /* Helper method. Creates random Cuboid based on objects parameters. Delegate shape creation to the standard
     * ShapeFactory */
    Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *randomShape();

public:
    ShapePair generate() override;
    std::string getDescription() const override;

    void setRadius(double _radius);
};

#endif // _BALL_FACTORY_H
