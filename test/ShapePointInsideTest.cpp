//--------------------------------------------------------------------------------------------
// Test of Cuboid::pointInside method
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "ShapePointInsideTest.h"
#include "utility/MockBC.h"
#include "../rsa3d/ShapeFactory.h"
#include "utility/BallFactory.h"
#include "../rsa3d/Utils.h"


namespace
{
    struct Results {
        std::size_t tested{};
        std::size_t overlapped{};
        std::size_t withPointInside{};
        std::size_t conflicts{};
        std::string factoryDesc;

        /* Prints test results onto given ostream */
        void print(std::ostream &ostr) {
            ostr << ">> Test results:" << std::endl;
            ostr << "factory           : " << this->factoryDesc << std::endl;
            ostr << "pairs tested      : " << this->tested << std::endl;
            ostr << "overlapped        : " << this->overlapped << std::endl;
            ostr << "with point inside : " << this->withPointInside << std::endl;
            ostr << "conflicts         : " << this->conflicts << std::endl;
        }
    };

    // Performs Cuboid::pointInside test. It generates a pair of cubes and checks results
    // given by Cuboid::overlap and Cuboid::pointInside. If pointInside gives true, so should
    // overlap do
    //----------------------------------------------------------------------------------------
    Results perform(ShapePairFactory &factory, unsigned long pairsToTest) {
        if (pairsToTest == 0) throw std::runtime_error("pairsToTest == 0");

        MockBC bc;
        Results results;
        results.tested = pairsToTest;
        results.factoryDesc = factory.getDescription();

        std::cout << ">> Starting..." << std::endl;
        for (std::size_t i = 0; i < pairsToTest; i++) {
            auto pair = factory.generate();
            bool overlap = pair.first()->overlap(&bc, pair.second());
            bool pi_first = pair.first()->pointInside(&bc, pair.second()->getPosition());
            bool pi_second = pair.second()->pointInside(&bc, pair.first()->getPosition());

            if (overlap)    results.overlapped++;
            if (pi_first)   results.withPointInside++;

            // pointInside and overlap conflict
            if (!overlap && (pi_first || pi_second))
                results.conflicts++;

            if ((i % 10000) == 9999)
                std::cout << (i + 1) << " pairs tested..." << std::endl;
        }
        return results;
    }
}

namespace shape_pitest
{
    int main(int argc, char **argv)
    {
        if (argc < 6)
            die("Usage: ./rsa shape_pitest [particle] [attr] [ball_radius] [max_tries]");

        double ballRadius = std::stod(argv[4]);
        unsigned long maxTries = std::stoul(argv[5]);
        if (ballRadius <= 0 || maxTries <= 0)
            die("Wrong input. Aborting.");

        ShapeFactory::initShapeClass(argv[2], argv[3]);
        BallFactory factory;
        factory.setRadius(ballRadius);

        Results results = perform(factory, maxTries);
        std::cout << std::endl;
        results.print(std::cout);

        return EXIT_SUCCESS;
    }
}