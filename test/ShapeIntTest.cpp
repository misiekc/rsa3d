//--------------------------------------------------------------------------------------------
// Test of cuboid intersection algorithm. It uses standard triangle-triangle intersection
// algorithm to check results given by Cuboid::overlap
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <typeinfo>
#include <memory>

#include "../rsa3d/Vector.h"
#include "../rsa3d/shapes/Cuboid.h"
#include "utility/ShapePairFactory.h"
#include "ShapeIntTest.h"
#include "../rsa3d/ShapeFactory.h"
#include "utility/BallFactory.h"
#include "../rsa3d/Utils.h"



namespace shape_inttest
{
    // Performs Cuboid::overlap algorithm check. It generates some random pairs of cuboids
    // and compares result given by ::MINE strategy and _strategy strategy
    //--------------------------------------------------------------------------------------------
    Results perform(ShapePairFactory *_factory, const OverlapStrategy *firstStrategy, const OverlapStrategy *secondStrategy,
                    int _max_tries)
    {
        shape_inttest::Results result;
        result.tries = _max_tries;
	
	    std::cout << ">> Performing " << typeid(*firstStrategy).name() << " and " << typeid(*secondStrategy).name()
	        << " for OverlapStrategy comparison..." << std::endl;
	    for (int i = 0; i < _max_tries; i++) {
	        ShapePairFactory::ShapePair pair = _factory->generate();

            bool firstIntersected = (bool)firstStrategy->overlap(pair.first, pair.second);
            bool secondIntersected = (bool)secondStrategy->overlap(pair.first, pair.second);

            if (firstIntersected)
                result.intersected++;
            else
                result.disjunct++;
	        
	        if (firstIntersected != secondIntersected) {  // Missed overlap return value, dump
	            result.missed++;
	            result.missed_dump.push_back(pair);
	        } else {  // Otherwisedelete pair
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
        if (argc < 6)
            die("Usage: ./rsa shape_inttest [particle] [attibutes] [ball_radius] [max_tries]");

        double ball_radius = std::stod(argv[4]);
        int max_tries = std::stoi(argv[5]);
        if (ball_radius <= 0 || max_tries <= 0)
            die("Wrong input. Aborting.");

        ShapeFactory::initShapeClass(argv[2], argv[3]);
        BallFactory * factory = BallFactory::getInstance();
        factory->setRadius(ball_radius);

        RND rnd;
        auto shape = ShapeFactory::createShape(&rnd);
        std::unique_ptr<const OverlapStrategyShape> osShape(dynamic_cast<const OverlapStrategyShape*>(shape));
        if (!osShape)
            die(std::string(argv[2]) + " is not OverlapStrategyShape");
        else if (osShape->getSupportedStrategies().size() < 2)
            die(std::string(argv[2]) + " does not supprot at least 2 strategies");

        auto strategies = osShape->getSupportedStrategies();
        std::unique_ptr<const OverlapStrategy> firstStrategy(osShape->createStrategy(strategies.front()));
        for (auto name = strategies.begin() + 1; name != strategies.end(); name++) {
            std::unique_ptr<const OverlapStrategy> secondStrategy(osShape->createStrategy(*name));
            Results results = perform(factory, firstStrategy.get(), secondStrategy.get(), max_tries);
            print_results(results);
            std::cout << std::endl;
            results.free_missed_pairs();
        }

        /*// Dump missed test to file
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
        }*/

        return EXIT_SUCCESS;
    }
}
