//--------------------------------------------------------------------------------------------
// Interface of factory producing pairs of random cuboids
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include "../shapes/Cuboid.h"


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
};
