//--------------------------------------------------------------------------------------------
// Interface of factory producing pairs of random cuboids
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CUBOID_PAIR_FACTORY_H
    #define _CUBOID_PAIR_FACTORY_H

#include <string>

#include "../../shapes/Cuboid.h"


class CuboidPairFactory
{
public:

    struct CuboidPair
    {
        Cuboid * first;
        Cuboid * second;
        
        CuboidPair()
        { }
        
        CuboidPair(Cuboid * first, Cuboid * second) : first(first), second(second)
        { }
        
        // Self-cleaning utility
        inline void free()
        {
            delete first;
            delete second;
        }
    };
    
    // Generates random cuboid pair for intersection checks
    virtual CuboidPair generate() = 0;
    virtual std::string getDescription() = 0;
};

#endif // _CUBOID_PAIR_FACTORY_H