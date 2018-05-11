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
#include "../rsa3d/ShapeFactory.h"
#include "../rsa3d/Utils.h"
#include "utility/BallFactory.h"
#include "utility/ShapePairFactory.h"

namespace
{
    /* A struct representing intersection test results */
    struct Results {
        unsigned long tries{};
        unsigned int intersected{};
        unsigned int disjunct{};
        std::vector<ShapePairFactory::ShapePair> missed_dump;

        std::size_t missed() const { return missed_dump.size(); }

        void report(ShapePairFactory::ShapePair &pair, bool firstIntersected, bool secondIntersected) {
            if (firstIntersected)   this->intersected++;
            else                    this->disjunct++;

            if (firstIntersected != secondIntersected)  // Missed overlap return value, dump
                this->missed_dump.push_back(pair);
        }

        /* Prints result info to given ostream */
        void print(std::ostream &ostr) const {
            ostr << ">> " << this->missed() << " from " << this->tries << " intersection results missed" << std::endl;
            ostr << ">> " << this->intersected << " shapes overlapped, " << this->disjunct;
            ostr << " shapes were disjunctive" << std::endl;
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
            if (results.missed() > 0) {
                if (!file.is_open())
                    file.open(filename);
                
                for (auto pair : results.missed_dump) {
                    file << "first = " << pair.first()->toWolfram() << ";" << std::endl;
                    file << "second = " << pair.second()->toWolfram() << ";" << std::endl;
                    file << "Graphics3D[{first, second}]" << std::endl;
                    file << std::endl;
                }
            }
        }
    };

    /* Creates a shape using shape factory and tries to cast it to OverlapStrategyShape */
    std::unique_ptr<RSAOverlapStrategyShape> acquire_shape() {
        RND rnd;
        auto shape = ShapeFactory::createShape(&rnd);
        auto osShape = dynamic_cast<RSAOverlapStrategyShape*>(shape);
        return std::unique_ptr<RSAOverlapStrategyShape>(osShape);
    }

    /* Check if a shape created by acquire_shape is valid for this test */
    void verifyShape(const std::string &name, const RSAOverlapStrategyShape *osShape) {
        if (osShape == nullptr)
            die(name + " is not OverlapStrategyShape");
        else if (osShape->getSupportedStrategies().size() < 2)
            die(name + " does not support at least 2 strategies");
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
    Results perform_test(ShapePairFactory &factory, const RSAOverlapStrategy &firstStrategy,
                         const RSAOverlapStrategy &secondStrategy, unsigned long maxTries) {
        Results result;
        result.tries = maxTries;

        std::cout << ">> Performing " << get_strategy_name(firstStrategy) << " and ";
        std::cout << get_strategy_name(secondStrategy) << " for OverlapStrategy comparison..." << std::endl;
        for (unsigned long i = 0; i < maxTries; i++) {
            ShapePairFactory::ShapePair pair = factory.generate();
            bool firstIntersected = (bool)firstStrategy.overlap(pair.first(), pair.second());
            bool secondIntersected = (bool)secondStrategy.overlap(pair.first(), pair.second());

            result.report(pair, firstIntersected, secondIntersected);
            if ((i % 10000) == 9999)
                std::cout << (i + 1) << " pairs tested..." << std::endl;
        }
        return result;
    }
}

namespace shape_ovtest
{
    int main(int argc, char **argv) {
        if (argc < 6)
            die("Usage: ./rsa_test shape_ovtest [particle] [attibutes] [ball_radius] [max_tries]");

        double ball_radius = std::stod(argv[4]);
        unsigned long max_tries = std::stoul(argv[5]);
        if (ball_radius <= 0 || max_tries <= 0)
            die("Wrong input. Aborting.");

        ShapeFactory::initShapeClass(argv[2], argv[3]);
        BallFactory factory;
        factory.setRadius(ball_radius);

        auto osShape = acquire_shape();
        verifyShape(argv[2], osShape.get());

        // Compare first strategy results with each following. Dump all missed pairs to a file
        auto strategies = osShape->getSupportedStrategies();
        auto firstStrategy = acquire_strategy(*osShape, strategies.front());
        MissedPairsDumper dumper("ovtest_dump.nb");
        for (auto name = strategies.begin() + 1; name != strategies.end(); name++) {
            auto secondStrategy = acquire_strategy(*osShape, *name);
            Results results = perform_test(factory, *firstStrategy, *secondStrategy, max_tries);
            results.print(std::cout);
            dumper.dump(results);
            if (results.missed() > 0) std::cout << ">> Missed pairs dumped to " << dumper.getFilename() << std::endl;
            std::cout << std::endl;
        }
        return EXIT_SUCCESS;
    }
}