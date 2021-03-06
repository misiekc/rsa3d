//--------------------------------------------------------------------------------------------
// Test of Cuboid::pointInside method
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "ShapePointInsideTest.h"
#include "../rsa3d/boundary_conditions/FreeBC.h"
#include "../rsa3d/shape/ShapeFactory.h"
#include "../rsa3d/shape/ConvexShape.h"
#include "utils/IndependentPairFactory.h"
#include "utils/InfoLooper.h"
#include "utils/UniformBallDistribution.h"
#include "utils/TestExitCodes.h"


namespace
{

	using ShapePair = RSAShapePairFactory::ShapePair;

	class Results {
	private:
	    RSAFreeBC bc;

	public:
        std::size_t tested{};
        std::size_t overlapped{};
        std::size_t withPointInside{};
        std::size_t conflicts{};
        std::string conflictExample;
        std::string factoryDesc;

        bool hasPointInside(const RSAShape *shape) {
            return dynamic_cast<const RSAConvexShape *>(shape) != nullptr;
        }

        auto isPointInside(const ShapePair &pair) {
            bool secondInsideFirst, firstInsideSecond;
            if (hasPointInside(pair.first())) {
                auto &firstConvex = dynamic_cast<const RSAConvexShape &>(*pair.first());
                auto &secondConvex = dynamic_cast<const RSAConvexShape &>(*pair.second());

                secondInsideFirst = firstConvex.pointInside(&bc, secondConvex.getPosition());
                firstInsideSecond = secondConvex.pointInside(&bc, firstConvex.getPosition());
            } else {
                RSAOrientation zeroOrientation;
                zeroOrientation.fill(0);

                secondInsideFirst = pair.first()->voxelInside(&bc, pair.second()->getPosition(), zeroOrientation, 0,
                        RSAShape::getAngularVoxelSize());
                firstInsideSecond = pair.second()->voxelInside(&bc, pair.first()->getPosition(), zeroOrientation, 0,
                        RSAShape::getAngularVoxelSize());
            }
            return std::make_pair(secondInsideFirst, firstInsideSecond);
        }

        void checkConflict(const ShapePair &pair, bool secondInsideFirst, bool firstInsideSecond, bool overlap) {
            if (!overlap && (secondInsideFirst || firstInsideSecond)) {
                if (conflicts == 0) {
                    std::stringstream sout;
                    pair.print(sout);
                    conflictExample = sout.str();
                }
                conflicts++;
            }
        }

        void evaluatePair(const ShapePair &pair) {
            auto [secondInsideFirst, firstInsideSecond] = isPointInside(pair);
            if (secondInsideFirst)
                withPointInside++;

            bool overlap = pair.first()->overlap(&bc, pair.second());
            if (overlap)
                overlapped++;

            checkConflict(pair, secondInsideFirst, firstInsideSecond, overlap);
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
        if (pairsToTest == 0)
            throw std::runtime_error("pairsToTest == 0");

        Results results;
        results.tested = pairsToTest;
        results.factoryDesc = factory.getDescription();

        std::cout << ">> Starting..." << std::endl;
        InfoLooper looper(pairsToTest, 10000, "pairs tested...");
        while(looper.step())
            results.evaluatePair(factory.generate());

        return results;
    }
}

namespace shape_pitest
{
    int main(int argc, char **argv)
    {
        if (argc < 6) {
            std::cerr << "Usage: ./rsa shape_pitest [particle] [attr] [ball_radius] [max_tries]" << std::endl;
            return TEST_ERROR;
        }

        double ballRadius = std::stod(argv[4]);
        unsigned long maxTries = std::stoul(argv[5]);
        if (ballRadius <= 0 || maxTries <= 0) {
            std::cerr << "Wrong input. Aborting." << std::endl;
            return TEST_ERROR;
        }

        ShapeFactory::initShapeClass(argv[2], argv[3]);
        RSAShape::setEarlyRejectionEnabled(false);
        UniformBallDistribution distribution(ballRadius);
        IndependentPairFactory factory(distribution);

        Results results = perform_test(factory, maxTries);
        std::cout << std::endl;
        results.print(std::cout);

        return results.conflicts == 0 ? TEST_SUCCESS : TEST_FAILURE;
    }
}
