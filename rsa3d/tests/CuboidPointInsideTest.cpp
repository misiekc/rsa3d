//--------------------------------------------------------------------------------------------
// Test of Cuboid::pointInside method
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "CuboidPointInsideTest.h"
#include "MockBC.h"

#include <iostream>

namespace cube_pitest
{   
    // Performs Cuboid::pointInside test. It generates a pair of cubes and checks results
    // given by Cuboid::overlap and Cuboid::pointInside. If pointInside gives true, so should
    // overlap do
    //----------------------------------------------------------------------------------------
    Results perform(CuboidPairFactory * _factory, std::size_t _pairs_to_test) {
        if (_pairs_to_test == 0)
            throw std::runtime_error("_pairs_to_test == 0");
        
        MockBC bc;
        Results results;
        results.tested = _pairs_to_test;
        results.factoryDesc = _factory->getDescription();
        
        std::cout << ">> Starting..." << std::endl;
        for (std::size_t i = 0; i < _pairs_to_test; i++) {
            CuboidPairFactory::CuboidPair pair;
	        pair = _factory->generate();
	        
	        bool overlap = pair.first->overlap(&bc, pair.second);
	        bool pi_first = pair.first->pointInside(&bc, pair.second->getPosition());
	        bool pi_second = pair.second->pointInside(&bc, pair.first->getPosition());
	        pair.free();
	        
	        if (overlap)
	            results.overlapped++;
	        if (pi_first)
	            results.withPointInside++;
	        
	        // pointInside and overlap conflict
	        if (!overlap && (pi_first || pi_second))
	            results.conflicts++;
	        
	        if ((i % 10000) == 9999)
	            std::cout << (i + 1) << " pairs tested..." << std::endl;
        }
        
        return results;
    }
    
    // Print pointInside test results
    //----------------------------------------------------------------------------------------
    void print_results(Results results) {
        std::cout << ">> Test results:" << std::endl;
        std::cout << "factory           : " << results.factoryDesc << std::endl;
        std::cout << "pairs tested      : " << results.tested << std::endl;
        std::cout << "overlapped        : " << results.overlapped << std::endl;
        std::cout << "with point inside : " << results.withPointInside << std::endl;
        std::cout << "conflicts         : " << results.conflicts << std::endl;
    }
}

void ble(){}
