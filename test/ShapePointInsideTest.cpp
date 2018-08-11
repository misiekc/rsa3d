//--------------------------------------------------------------------------------------------
// Test of Cuboid::pointInside method
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "ShapePointInsideTest.h"
#include "utility/MockBC.h"
#include "../rsa3d/shape/ShapeFactory.h"
#include "utility/IndependentPairFactory.h"
#include "../rsa3d/Utils.h"
#include "utility/InfoLooper.h"
#include "../rsa3d/shape/ConvexShape.h"
#include "utility/UniformBallDistribution.h"


namespace
{

	using ShapePair = RSAShapePairFactory::ShapePair;

	struct Results {
        std::size_t tested{};
        std::size_t overlapped{};
        std::size_t withPointInside{};
        std::size_t conflicts{};
        std::string conflictExample;
        std::string factoryDesc;

        /* Prints Graphics3D with shapes pair on the out ostream */
        void printPair(const ShapePair &pair, std::ostream &out) const {
            out << "Graphics3D[{ " << std::endl;
            out << pair.first()->toWolfram() << ", " << std::endl;
            out << pair.second()->toWolfram() << std::endl;
            out << "}]";
        }

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

            if (this->conflicts > 0) {
            	ostr << "[INFO] Example conflict configuration: two non-overlapping ";
            	ostr << "(according to the overlap() method) shapes." << std::endl;
            	ostr << "[INFO] One of them has a center inside the second pointInside() zone:";
            	ostr << std::endl;
            	ostr << this->conflictExample;
            	ostr << std::endl;
            }
        }
    };

    Results perform_test(RSAShapePairFactory &factory, unsigned long pairsToTest) {
        if (pairsToTest == 0) throw std::runtime_error("pairsToTest == 0");

        RSAMockBC bc;
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
            if (!overlap && (pi_first || pi_second)){
            	if (results.conflicts == 0) {
            		std::stringstream sout;
            		results.printPair(pair, sout);
            		results.conflictExample = sout.str();
            	}
                results.conflicts++;
            }
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
        UniformBallDistribution distribution(ballRadius);
        IndependentPairFactory factory(distribution);

        Results results = perform_test(factory, maxTries);
        std::cout << std::endl;
        results.print(std::cout);

        return EXIT_SUCCESS;
    }
}
