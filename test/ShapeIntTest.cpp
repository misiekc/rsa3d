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

namespace
{
    /* A struct representing intersection test results */
    struct Results {
        int tries = 0;
        int missed = 0;
        int intersected = 0;
        int disjunct = 0;
        std::vector<ShapePairFactory::ShapePair> missed_dump;

        inline void free_missed_pairs() {
            for (auto pair : missed_dump)
                pair.free();
        }
    };

    /* Creates a shape using shape factory and tries to cast it to OverlapStrategyShape */
    std::unique_ptr<RSAOverlapStrategyShape> acquireShape() {
        RND rnd;
        auto shape = ShapeFactory::createShape(&rnd);
        auto osShape = dynamic_cast<RSAOverlapStrategyShape*>(shape);
        return std::unique_ptr<RSAOverlapStrategyShape>(osShape);
    }

    /* Check if a shape created by aquireShape is valid fpor this test */
    void verifyShape(const std::string &name, const RSAOverlapStrategyShape *osShape) {
        if (osShape == nullptr)
            die(name + " is not OverlapStrategyShape");
        else if (osShape->getSupportedStrategies().size() < 2)
            die(name + " does not support at least 2 strategies");
    }

    /* Performs a comparison of two strategies */
    Results perform(ShapePairFactory &factory, const RSAOverlapStrategy &firstStrategy,
                    const RSAOverlapStrategy &secondStrategy, int _max_tries)
    {
        Results result;
        result.tries = _max_tries;

        std::cout << ">> Performing " << typeid(firstStrategy).name() << " and " << typeid(secondStrategy).name()
                  << " for OverlapStrategy comparison..." << std::endl;
        for (int i = 0; i < _max_tries; i++) {
            ShapePairFactory::ShapePair pair = factory.generate();

            bool firstIntersected = (bool)firstStrategy.overlap(pair.first, pair.second);
            bool secondIntersected = (bool)secondStrategy.overlap(pair.first, pair.second);

            if (firstIntersected)
                result.intersected++;
            else
                result.disjunct++;

            if (firstIntersected != secondIntersected) {  // Missed overlap return value, dump
                result.missed++;
                result.missed_dump.push_back(pair);
            } else {  // Otherwise delete pair
                pair.free();
            }

            if ((i % 10000) == 9999)
                std::cout << (i + 1) << " pairs tested..." << std::endl;
        }

        return result;
    }

    /* Prints result to a standard output */
    void print_results(Results results)
    {
        std::cout << ">> " << results.missed << " from " << results.tries << " intersection results missed" << std::endl;
        std::cout << ">> " << results.intersected << " shapes overlapped, " << results.disjunct << " shapes were disjunctive" << std::endl;
    }

    /* Dumps missed pair onto std::ostream ostr */
    void dump_missed_pairs(Results results, std::ostream &ostr)
    {
        for (auto pair : results.missed_dump) {
            ostr << "first = " << pair.first->toWolfram() << ";" << std::endl;
            ostr << "second = " << pair.second->toWolfram() << ";" << std::endl;
            ostr << "Graphics3D[{first, second}]" << std::endl;
            ostr << std::endl;
        }
    }

    /* Dumps missed pairs to a file of given name */
    void dump_missed_pairs_to_file(Results results, std::ostream &dumpFile) {
        if (results.missed > 0) {
            dump_missed_pairs(results, dumpFile);
            std::cout << ">> Missed pairs dumped to file " << std::endl;
        }
    }
}

namespace shape_inttest
{
    int main(int argc, char **argv)
    {
        if (argc < 6)
            die("Usage: ./rsa_test shape_inttest [particle] [attibutes] [ball_radius] [max_tries]");

        double ball_radius = std::stod(argv[4]);
        int max_tries = std::stoi(argv[5]);
        if (ball_radius <= 0 || max_tries <= 0)
            die("Wrong input. Aborting.");

        ShapeFactory::initShapeClass(argv[2], argv[3]);
        BallFactory factory;
        factory.setRadius(ball_radius);

        std::ofstream dumpFile("inttest_dump.nb");
        if (!dumpFile)  die("Cannot open inttest_dump.nb to write");

        auto osShape = acquireShape();
        verifyShape(argv[2], osShape.get());

        // Compare first strategy results with each following. Dump all missed pairs to a file
        auto strategies = osShape->getSupportedStrategies();
        std::unique_ptr<const RSAOverlapStrategy> firstStrategy(osShape->createStrategy(strategies.front()));
        for (auto name = strategies.begin() + 1; name != strategies.end(); name++) {
            std::unique_ptr<const RSAOverlapStrategy> secondStrategy(osShape->createStrategy(*name));
            Results results = perform(factory, *firstStrategy, *secondStrategy, max_tries);
            print_results(results);
            dump_missed_pairs_to_file(results, dumpFile);
            std::cout << std::endl;
            results.free_missed_pairs();
        }
        dumpFile.close();
        return EXIT_SUCCESS;
    }
}
