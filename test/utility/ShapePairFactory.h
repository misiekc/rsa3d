//--------------------------------------------------------------------------------------------
// Interface of factory producing pairs of random cuboids
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CUBOID_PAIR_FACTORY_H
    #define _CUBOID_PAIR_FACTORY_H

#include <string>

#include "../../rsa3d/shapes/Cuboid.h"
#include "../../rsa3d/Parameters.h"


class ShapePairFactory
{
public:

    struct ShapePair
    {
        Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> * first;
        Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> * second;
        
        ShapePair()
        { }
        
        ShapePair(Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> * first, Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> * second) : first(first), second(second)
        { }
        
        // Self-cleaning utility
        inline void free()
        {
            delete first;
            delete second;
        }
    };
    
    // Generates random cuboid pair for intersection checks
    virtual ShapePair generate() = 0;
    virtual std::string getDescription() = 0;
};

#endif // _CUBOID_PAIR_FACTORY_H
