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
#include <vector>
#include <fstream>

#include "../Vector.h"
#include "../Intersection.h"
#include "../shapes/Cuboid.h"
#include "../ShapeFactory.h"
#include "MockBC.h"


// Helper methods
//--------------------------------------------------------------------------------------------
namespace
{   
	int             no = 0;
    
    // Struct representing intersection test results
    //--------------------------------------------------------------------------------------------
    struct TestResult {
        int tries = 0;
        int missed = 0;
        int intersected = 0;
        int disjunct = 0;
    };
    
    
    // Helper method. Creates random Cuboid based on global parameters
    //--------------------------------------------------------------------------------------------
    Cuboid * random_cuboid(RND * rnd, double box_halfsize)
    {
	    double trans[3];
	    Cuboid * cube = (Cuboid*)ShapeFactory::createShape(rnd);
	    cube->no = no++;
	    trans[0] = (rnd->nextValue() * 2 - 1) * box_halfsize;
	    trans[1] = (rnd->nextValue() * 2 - 1) * box_halfsize;
	    trans[2] = (rnd->nextValue() * 2 - 1) * box_halfsize;
	    cube->translate(trans);
	    return cube;
    } 
    
    // Helper method. Prints intersection test results onto the standard output
    //--------------------------------------------------------------------------------------------
    void print_test_result(TestResult result)
    {
	    std::cout << ">> " << result.missed << " from " << result.tries << " intersection results missed" << std::endl;
	    std::cout << ">> " << result.intersected << " cuboids overlapped, " << result.disjunct << " cuboids were disjunctive" << std::endl;
    }
    
    // Helper method. Performs comparison of Cuboid::OverlapStrategy::MINE and the second from
    // parameter
    //--------------------------------------------------------------------------------------------
    TestResult perform_strategy_comparison(Cuboid::OverlapStrategy _second, double box_halfsize, int max_tries)
    {
        Cuboid  *cube1, *cube2;
        TestResult result;
	    bool    mine_intersected, second_intersected;
	    MockBC  bc;
        RND     rnd;
        std::stringstream missed_dump_stream;
	
        result.tries = max_tries;
	
	    for (int i = 0; i < max_tries; i++) {
	        cube1 = random_cuboid(&rnd, box_halfsize);
	        cube2 = random_cuboid(&rnd, box_halfsize);
	        
	        Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy::MINE);
	        mine_intersected = (bool)cube1->overlap(&bc, cube2);
	        Cuboid::setOverlapStrategy(_second);
	        second_intersected = (bool)cube1->overlap(&bc, cube2);
	        
	        if (mine_intersected != second_intersected) {
	            result.missed++;
	            missed_dump_stream << cube1->toWolfram() << std::endl;
	            missed_dump_stream << cube2->toWolfram() << std::endl;
	            missed_dump_stream << "pair" << result.missed << " = {cube" << cube1->no << ", cube" << cube2->no << "};" << std::endl;  
	            missed_dump_stream << "Graphics3D[pair" << result.missed << "]";
	            missed_dump_stream << std::endl;
	        } else {
	            if (mine_intersected)
	                result.intersected++;
	            else
	                result.disjunct++;
	        }
	        
	        if ((i % 10000) == 9999)
	            std::cout << (i + 1) << " pairs tested..." << std::endl;
	        
	        delete cube1;
	        delete cube2;
	    }
	    
	    // Dump missed test to stdout and file
	    if (result.missed > 0) {
	        if (result.missed < 10)
	            std::cout << missed_dump_stream.str();
	        std::ofstream dump_file("inttest_dump.nb");
	        if (dump_file) {
	            dump_file << missed_dump_stream.str();
	            dump_file.close();
	            std::cout << ">> Missed pairs dumped to inttest_dump.txt" << std::endl;
	        } else {
	            std::cout << ">> Could not write to inttest_dump.txt";
	        }
	    }
	    return result;
    }
} 


namespace cube_inttest
{

    // Runs intersection::tri_tri3D selftest. Two simple almoast-intersecting and 
    // almoast-not-intersecting cases are checked
    //--------------------------------------------------------------------------------------------
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
    // and compares result given by different overlap strategies
    //--------------------------------------------------------------------------------------------
    void perform(double sizex, double sizey, double sizez, double box_halfsize, int max_tries)
    {
        // Test intersection::tri_tri3D
        TriTriInt_selftest_run();       
        
        std::stringstream stream;
        stream << "3 " << sizex << " " <<  sizey << " " << sizez;
	    ShapeFactory::initShapeClass("Cuboid", stream.str());
	
	    TestResult result;
	
	    /*std::cout << ">> Performing ::MINE and ::TRI_TRI for Cuboid::OverlapStrategy comparison..." << std::endl;
	    result = perform_strategy_comparison(Cuboid::OverlapStrategy::TRI_TRI, box_halfsize, max_tries);
	    print_test_result(result);
	    
	    std::cout << std::endl;*/
	    
	    std::cout << ">> Performing ::MINE and ::SAT for Cuboid::OverlapStrategy comparison..." << std::endl;
        result = perform_strategy_comparison(Cuboid::OverlapStrategy::SAT, box_halfsize, max_tries);
	    print_test_result(result);
    }
}
