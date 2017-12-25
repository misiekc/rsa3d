//--------------------------------------------------------------------------------------------
// Test of Cuboid::pointInside method
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CUBOID_POINT_INSDE_TEST_H
    #define _CUBOID_POINT_INSDE_TEST_H
    

#include "utility/CuboidPairFactory.h"

#include <string>


namespace cube_pitest
{
    struct Results {
        std::size_t tested = 0;
        std::size_t overlapped = 0;
        std::size_t withPointInside = 0;
        std::size_t conflicts = 0;
        std::string factoryDesc;
    };
    
    Results perform(CuboidPairFactory * _factory, std::size_t _pairs_to_test);
    void print_results(Results results);
    void main(int argc, char **argv);
}

void ble();
#endif // _CUBOID_POINT_INSDE_TEST_H
