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
#include <memory>
#include <utility>
#include <typeinfo>
#include <cxxabi.h>

#include "ShapeOverlapTest.h"
#include "../rsa3d/shape/ShapeFactory.h"
#include "../rsa3d/utils/Utils.h"
#include "utils/IndependentPairFactory.h"
#include "utils/ShapePairFactory.h"
#include "utils/InfoLooper.h"
#include "utils/ParallelPairFactory.h"
#include "utils/UniformBallDistribution.h"
#include "utils/TestExitCodes.h"

namespace
{
    /* A struct representing intersection test results */
    struct Results {
        unsigned long tries{};
        unsigned int intersected{};
        unsigned int disjunct{};
        std::vector<RSAShapePairFactory::ShapePair> missed_dump;

        std::size_t missed() const { return missed_dump.size(); }

        void report(RSAShapePairFactory::ShapePair &pair, bool firstIntersected, bool secondIntersected) {
            if (firstIntersected)   this->intersected++;
            else                    this->disjunct++;

            if (firstIntersected != secondIntersected)  // Missed overlap return value, dump
                this->missed_dump.push_back(pair);
        }

        /* Prints result info to given ostream */
        void print(std::ostream &ostr) const {
            ostr << "[INFO] " << this->intersected << " shapes overlapped, " << this->disjunct;
            ostr << " shapes were disjunctive" << std::endl;
            ostr << (this->missed() == 0 ? "[SUCCESS] " : "[FAILURE] ");
            ostr << this->missed() << " from " << this->tries << " intersection results missed" << std::endl;
        }
    };

    /* Helper class for dumping missed pairs. Creates dump file only if it is needed. */
    class MissedPairsDumper {
    private:
        std::string filename;
        std::ofstream file;

    public:
        explicit MissedPairsDumper(std::string filename) : filename(std::move(filename)) {}
        
        ~MissedPairsDumper() { 
            if (file.is_open()) file.close();
        }

        const std::string &getFilename() const { return filename; }

        void dump(Results results) {
            if (results.missed() <= 0)  return;

            if (!file.is_open())
                file.open(filename);
            for (auto pair : results.missed_dump) {
                file << "first = " << pair.first()->toWolfram() << ";" << std::endl;
                file << "second = " << pair.second()->toWolfram() << ";" << std::endl;
                file << "Graphics3D[{first, second}]" << std::endl;
                file << std::endl;
            }
        }
    };

    /* Creates a shape using shape factory and tries to cast it to OverlapStrategyShape */
    std::unique_ptr<RSAOverlapStrategyShape> acquire_overlap_strategy_shape() {
        RND rnd;
        auto shape = ShapeFactory::createShape(&rnd);
        auto osShape = dynamic_cast<RSAOverlapStrategyShape*>(shape);
        return std::unique_ptr<RSAOverlapStrategyShape>(osShape);
    }

    /* Check if a shape created by acquire_overlap_strategy_shape is valid for this test */
    bool is_shape_supported(const std::string &name, const RSAOverlapStrategyShape *osShape) {
        if (osShape == nullptr) {
            std::cerr << name << " is not OverlapStrategyShape" << std::endl;
            return false;
        } else if (osShape->getSupportedStrategies().size() < 2) {
            std::cerr << name << " does not support at least 2 strategies" << std::endl;
            return false;
        } else {
            return true;
        }
    }

    /* Wraps some ugliness */
    std::unique_ptr<const RSAOverlapStrategy> acquire_strategy(const RSAOverlapStrategyShape &shape,
                                                               const std::string &strategyName) {
        return std::unique_ptr<const RSAOverlapStrategy>(shape.createStrategy(strategyName));
    }

    /* Returns formatted (demangled) strategy name based on typeid and abi's demangling */
    std::string get_strategy_name(const RSAOverlapStrategy &strategy) {
        int status;
        auto _name = abi::__cxa_demangle(typeid(strategy).name(), nullptr, nullptr, &status);
        if (status != 0)    throw std::runtime_error("demangle error");
        std::string name(_name);
        std::free(_name);
        return name;
    }

    /* Performs a comparison of two strategies */
    Results test_strategy_pair(RSAShapePairFactory &factory, const RSAOverlapStrategy &firstStrategy,
                               const RSAOverlapStrategy &secondStrategy, unsigned long maxTries) {
        Results result;
        result.tries = maxTries;

        std::cout << "[INFO] Performing " << get_strategy_name(firstStrategy) << " and ";
        std::cout << get_strategy_name(secondStrategy) << " for OverlapStrategy comparison..." << std::endl;
        InfoLooper looper(maxTries, 25000, "pairs tested...");
        while (looper.step()) {
            RSAShapePairFactory::ShapePair pair = factory.generate();
            bool firstIntersected = firstStrategy.overlap(pair.first(), pair.second());
            bool secondIntersected = secondStrategy.overlap(pair.first(), pair.second());
            result.report(pair, firstIntersected, secondIntersected);
        }
        return result;
    }

    /* Performs test of all strategies, prints info and dumps missed pairs. Returns true on success. */
    bool perform_test(const RSAOverlapStrategyShape &osShape, RSAShapePairFactory &factory, unsigned long maxTries,
                      const std::string &dumpFile) {
        auto strategies = osShape.getSupportedStrategies();
        auto firstStrategy = acquire_strategy(osShape, strategies.front());
        MissedPairsDumper dumper(dumpFile);
        Results results;
        for (auto name = strategies.begin() + 1; name != strategies.end(); name++) {
            auto secondStrategy = acquire_strategy(osShape, *name);
            results = test_strategy_pair(factory, *firstStrategy, *secondStrategy, maxTries);
            results.print(std::cout);
            dumper.dump(results);
            if (results.missed() > 0) std::cout << "[INFO] Missed pairs dumped to " << dumper.getFilename() << std::endl;
            std::cout << std::endl;
        }
        return (results.missed() == 0);
    }
}

namespace shape_ovtest
{
    int main(int argc, char **argv) {
        if (argc < 6) {
            std::cerr << "Usage: ./rsa_test shape_ovtest [particle] [attibutes] [ball_radius] [max_tries]" << std::endl;
            return TEST_ERROR;
        }

        double ball_radius = std::stod(argv[4]);
        unsigned long max_tries = std::stoul(argv[5]);
        if (ball_radius <= 0 || max_tries <= 0) {
            std::cerr << "Wrong input. Aborting." << std::endl;
            return TEST_ERROR;
        }

        ShapeFactory::initShapeClass(argv[2], argv[3]);
        RSAShape::setEarlyRejectionEnabled(false);

        auto osShape = acquire_overlap_strategy_shape();
        if (!is_shape_supported(argv[2], osShape.get()))
            return TEST_UNSUPPORTED;

        UniformBallDistribution distribution(ball_radius);
        IndependentPairFactory factory(distribution);
        std::cout << "[INFO] Performing test with unoriented shapes -------------------------------------" << std::endl;
        if (!perform_test(*osShape, factory, max_tries, "ovtest_anisotropic_dump.nb"))
            return TEST_FAILURE;

        ParallelPairFactory isotropicFactory(distribution);
        std::cout << "[INFO] Performing test with oriented shapes ---------------------------------------" << std::endl;
        if (!perform_test(*osShape, isotropicFactory, max_tries, "ovtest_isotropic_dump.nb"))
            return TEST_FAILURE;

        return TEST_SUCCESS;
    }
}