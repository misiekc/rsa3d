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
#include "utility/InfoLooper.h"
#include "../rsa3d/ConvexShape.h"


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
            std::size_t percentAccuracy = 100 * this->withPointInside / this->overlapped;

            ostr << ">> Test results:" << std::endl;
            ostr << "factory               : " << this->factoryDesc << std::endl;
            ostr << "pairs tested          : " << this->tested << std::endl;
            ostr << "overlapped            : " << this->overlapped << std::endl;
            ostr << "with point inside     : " << this->withPointInside << std::endl;
            ostr << "point inside accuracy : " << percentAccuracy << "%" << std::endl;
            ostr << "conflicts             : " << this->conflicts << std::endl;
        }
    };

    Results perform_test(RSAShapePairFactory &factory, unsigned long pairsToTest) {
        if (pairsToTest == 0) throw std::runtime_error("pairsToTest == 0");

        MockBC bc;
        Results results;
        results.tested = pairsToTest;
        results.factoryDesc = factory.getDescription();

        std::cout << ">> Starting..." << std::endl;
        InfoLooper looper(pairsToTest, 10000, "pairs tested...");
        while(looper.step()) {
            auto pair = factory.generate();
            auto &first = dynamic_cast<RSAConvexShape&>(*pair.first());
            auto &second = dynamic_cast<RSAConvexShape&>(*pair.second());
            
            bool overlap = first.overlap(&bc, &second);
            bool pi_first = first.pointInside(&bc, second.getPosition());
            bool pi_second = second.pointInside(&bc, first.getPosition());

            if (overlap)    results.overlapped++;
            if (pi_first)   results.withPointInside++;

            // pointInside and overlap conflict
            if (!overlap && (pi_first || pi_second))
                results.conflicts++;
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

        Results results = perform_test(factory, maxTries);
        std::cout << std::endl;
        results.print(std::cout);

        return EXIT_SUCCESS;
    }
}