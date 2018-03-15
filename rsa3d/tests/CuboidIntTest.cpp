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
#include <fstream>

#include "../Vector.h"
#include "../Intersection.h"
#include "../shapes/Cuboid.h"
#include "utility/MockBC.h"
#include "utility/ShapePairFactory.h"
#include "CuboidIntTest.h"
#include "../shapes/cube_strategies/MineOverlap.h"
#include "../ShapeFactory.h"
#include "utility/BallFactory.h"
#include "../shapes/cube_strategies/SATOverlap.h"
#include "../shapes/cube_strategies/TriTriOverlap.h"
#include "../shapes/cube_strategies/OptimizedSATOverlap.h"
#include "../Utils.h"


namespace
{
    MineOverlap mineOverlap;
}

namespace cube_inttest
{
    // Performs Cuboid::overlap algorithm check. It generates some random pairs of cuboids
    // and compares result given by ::MINE strategy and _strategy strategy
    //--------------------------------------------------------------------------------------------
    Results perform(ShapePairFactory * _factory, OverlapStrategy * _strategy, int _max_tries)
    {
        cube_inttest::Results result;
	    bool    mine_intersected, second_intersected;
	    MockBC  bc;
        RND     rnd;
	
        result.tries = _max_tries;
	
	    std::cout << ">> Performing " << mineOverlap.getName() << " and " << _strategy->getName()
	        << " for OverlapStrategy comparison..." << std::endl;
	    for (int i = 0; i < _max_tries; i++) {
	        ShapePairFactory::ShapePair pair = _factory->generate();
	        
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
            _ostr << ((Cuboid*)(pair.first))->toWolfram() << std::endl;
            _ostr << ((Cuboid*)(pair.second))->toWolfram() << std::endl;
            _ostr << "pair" << pair_no << " = {cube" << pair.first->no << ", cube" << pair.second->no << "};" << std::endl;  
            _ostr << "Graphics3D[pair" << pair_no << "]" << std::endl;
            _ostr << std::endl;
            pair_no++;
	    }
    }

    // Performs cuboid intersection test with parameters passed to process
    //----------------------------------------------------------------------------------------
    int main(int argc, char **argv)
    {
        if (argc < 7)
            die("Usage: ./rsa cube_inttest [size_x] [size_y] [size_z] [ball_radius] [max_tries]");

        double ball_radius = std::stod(argv[5]);
        int max_tries = std::stoi(argv[6]);
        if (ball_radius <= 0 || max_tries <= 0)
            die("Wrong input. Aborting.");

        std::stringstream init_stream;
        init_stream << "3 " << argv[2] << " " << argv[3] << " " << argv[4];
        ShapeFactory::initShapeClass("Cuboid", init_stream.str());
        BallFactory * factory = BallFactory::getInstance();
        factory->setRadius(ball_radius);

        // Test SAT
        SATOverlap satOverlap;
        cube_inttest::Results results = cube_inttest::perform(factory, &satOverlap, max_tries);
        print_results(results);
        std::cout << std::endl;
        results.free_missed_pairs();

        // Test TRI_TRI
        TriTriOverlap triTriOverlap;
        results = cube_inttest::perform(factory, &triTriOverlap, max_tries);
        print_results(results);
        std::cout << std::endl;
        results.free_missed_pairs();

        // Test optimised SAT
        OptimizedSATOverlap optimizedSATOverlap;
        results = cube_inttest::perform(factory, &optimizedSATOverlap, max_tries);
        print_results(results);
        std::cout << std::endl;

        // Dump missed test to file
        if (results.missed > 0) {
            std::ofstream dump_file("inttest_dump.nb");
            if (dump_file) {
                dump_missed_pairs(results, dump_file);
                results.free_missed_pairs();
                dump_file.close();
                std::cout << ">> Missed pairs dumped to inttest_dump.nb" << std::endl;
            } else {
                std::cout << ">> Could not write to inttest_dump.nb";
            }
        }

        return EXIT_SUCCESS;
    }
}
