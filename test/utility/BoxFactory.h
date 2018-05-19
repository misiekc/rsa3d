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
    std::array<double, RSA_SPATIAL_DIMENSION> halfsize{};
    unsigned int no{};

    /* Helper method. Creates random Cuboid based on objects parameters. Delegate cuboid creation to standard
     * ShapeFactory */
    RSAShape *randomShape();

public:
    BoxFactory() : BoxFactory(halfsize) {};
    explicit BoxFactory(const std::array<double, RSA_SPATIAL_DIMENSION> &halfsize) : halfsize(halfsize) {};
    explicit BoxFactory(double halfsize);

    std::string getDescription() const override;
    ShapePair generate() override;

    void setBoxSize(const std::array<double, RSA_SPATIAL_DIMENSION> &halfsize);
};

#endif // _BOX_FACTORY_H
