//--------------------------------------------------------------------------------------------
// Test of cuboid intersection algorithm. It uses standard triangle-triangle intersection
// algorithm to check results given by Cuboid::overlap
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include <iostream>
#include <cassert>
#include <cstring>
#include <sstream>

#include "../Vector.h"
#include "../Intersection.h"
#include "../shapes/Cuboid.h"
#include "../ShapeFactory.h"
#include "MockBC.h"
#include "CuboidPairFactory.h"
#include "CuboidIntTest.h" 

namespace
{
    MineOverlap mineOverlap;
}

namespace cube_inttest
{

    // Runs intersection::tri_tri3D selftest. Two simple almoast-intersecting and 
    // almoast-not-intersecting cases are checked
    //----------------------------------------------------------------------------------------
    void TriTriInt_selftest_run()
    {
        Vector<3> triangle1[] = {
            Vector<3>{{0, -1, 0}},
            Vector<3>{{0, 1, 0}},
            Vector<3>{{0, 0, 1}}
        };
        
        // Triangle on z = 1.00000001 plane containing x = 0, y = 0
        Vector<3> triangle2[] = {
            Vector<3>{{1.4, 0.8, 1.00000001}},
            Vector<3>{{-1.3, 0, 1.00000001}},
            Vector<3>{{0, -4.5, 1.00000001}}
        };
        
        // Triangle on z = 0.99999999 plane containing x = 0, y = 0
        Vector<3> triangle3[] = {
            Vector<3>{{1.4, 0.8, 0.99999999}},
            Vector<3>{{-1.3, 0, 0.99999999}},
            Vector<3>{{0, -4.5, 0.99999999}}
        };
        
        Vector<3> triangle4[] = {
            Vector<3>{{-1, -1, 0}},
            Vector<3>{{1, 0, 2.00000001}},
            Vector<3>{{-1, 1, 0}}
        };
        
        Vector<3> triangle5[] = {
            Vector<3>{{-1, -1, 0}},
            Vector<3>{{1, 0, 1.99999999}},
            Vector<3>{{-1, 1, 0}}
        };
        
        std::cout << std::boolalpha;
        std::cout << ">> Performing intersection::tri_tri3D quick selftest..." << std::endl;
        assert(intersection::tri_tri3D(triangle1, triangle2) == false);
        assert(intersection::tri_tri3D(triangle1, triangle3) == true);
        assert(intersection::tri_tri3D(triangle1, triangle4) == false);
        assert(intersection::tri_tri3D(triangle1, triangle5) == true);
        assert(intersection::tri_tri3D(triangle2, triangle1) == false);
        assert(intersection::tri_tri3D(triangle3, triangle1) == true);
        assert(intersection::tri_tri3D(triangle4, triangle1) == false);
        assert(intersection::tri_tri3D(triangle5, triangle1) == true);
        std::cout << "Passed." << std::endl;
    }


    // Performs Cuboid::overlap algorithm check. It generates some random pairs of cuboids
    // and compares result given by ::MINE strategy and _strategy strategy
    //--------------------------------------------------------------------------------------------
    Results perform(CuboidPairFactory * _factory, OverlapStrategy * _strategy, int _max_tries)
    {
        cube_inttest::Results result;
	    bool    mine_intersected, second_intersected;
	    MockBC  bc;
        RND     rnd;
	
        result.tries = _max_tries;
	
	    std::cout << ">> Performing " << mineOverlap.getName() << " and " << _strategy->getName()
	        << " for OverlapStrategy comparison..." << std::endl;
	    for (int i = 0; i < _max_tries; i++) {
	        CuboidPairFactory::CuboidPair pair = _factory->generate();
	        
	        Cuboid::setOverlapStrategy(&mineOverlap);
	        mine_intersected = (bool)pair.first->overlap(&bc, pair.second);
	        Cuboid::setOverlapStrategy(_strategy);
	        second_intersected = (bool)pair.first->overlap(&bc, pair.second);
	        
	        if (mine_intersected != second_intersected) {  // Missed overlap return value, dump
	            result.missed++;
	            result.missed_dump.push_back(pair);
	        } else {  // Otherwise, update statistics and delete pair
	            if (mine_intersected)
	                result.intersected++;
	            else
	                result.disjunct++;
	            pair.free();
	        }
	        
	        if ((i % 10000) == 9999)
	            std::cout << (i + 1) << " pairs tested..." << std::endl;  
	    }
	    
	    return result;
    }
    
    // Prints intersection test results onto the standard output
    //--------------------------------------------------------------------------------------------
    void print_results(Results _results)
    {
	    std::cout << ">> " << _results.missed << " from " << _results.tries << " intersection results missed" << std::endl;
	    std::cout << ">> " << _results.intersected << " cuboids overlapped, " << _results.disjunct << " cuboids were disjunctive" << std::endl;
    }
    
    // Dumps missed pair onto the stream in the form of Wolfram code to display 3d view
    //--------------------------------------------------------------------------------------------
    void dump_missed_pairs(Results _results, std::ostream & _ostr)
    {
        std::size_t pair_no = 0;
        for (auto pair : _results.missed_dump) {
            _ostr << pair.first->toWolfram() << std::endl;
            _ostr << pair.second->toWolfram() << std::endl;
            _ostr << "pair" << pair_no << " = {cube" << pair.first->no << ", cube" << pair.second->no << "};" << std::endl;  
            _ostr << "Graphics3D[pair" << pair_no << "]" << std::endl;
            _ostr << std::endl;
            pair_no++;
	    }
    }
}
