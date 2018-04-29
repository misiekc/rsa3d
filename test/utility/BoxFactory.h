//--------------------------------------------------------------------------------------------
// Factory generating cuboids from box of given size
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _BOX_FACTORY_H
    #define _BOX_FACTORY_H

#include "ShapePairFactory.h"


/**
 * @brief ShapePairFactory generating shapes from a box of given size.
 */
class BoxFactory : public ShapePairFactory
{
private:
    double halfsizeX{0.5};
    double halfsizeY{0.5};
    double halfsizeZ{0.5};
    unsigned int no{};
    RND rnd{};

    /* Helper method. Creates random Cuboid based on objects parameters. Delegate cuboid creation to standard
     * ShapeFactory */
    Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *randomShape();

public:
    std::string getDescription() const override;
    ShapePair generate() override;

    void setBoxSize(double _halfsize_x, double _halfsize_y, double _halfsize_z);
};

#endif // _BOX_FACTORY_H
