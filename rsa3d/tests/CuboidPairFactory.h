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
    };
    
    // Generates random cuboid pair for intersection checks
    virtual CuboidPair generate() = 0;
};
