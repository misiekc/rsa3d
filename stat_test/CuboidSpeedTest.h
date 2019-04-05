//--------------------------------------------------------------------------------------------
// Comparison of different overlap algorithms speed depending on overlap probability
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CUBOID_SPEED_TEST_H
    #define _CUBOID_SPEED_TEST_H


#include "../rsa3d/Config.h"
#include "../rsa3d/shape/ShapeFactory.h"
#include "utility/IndependentPairFactory.h"
#include "utility/ShapePairFactory.h"
#include "../rsa3d/Quantity.h"
#include "../rsa3d/shape/shapes/cuboid/CuboidOverlapStrategy.h"

#include <vector>
#include <ostream>


namespace cube_speedtest
{
    /* Test context */
    struct Context {
        std::size_t pairs{};
        std::size_t repeats{};
        std::vector<CuboidOverlapStrategy *> strategies;
        std::vector<double> ballRadia;

        void load(const std::string &_filename);
        ~Context();
    };

    /* End results for single strategy */
    struct StrategyResult {
        CuboidOverlapStrategy * strategy;
        Quantity        time;
    };

    /* End results for single pair factory size */
    struct Result {
        Quantity        numOverlapped;
        std::size_t     numAll{};
        Quantity        overlapProb;
        std::vector<StrategyResult>   strategyResults;
        std::string     factoryDesc;

        void print() const;
        static void printVector(const std::vector<Result> &_results);
        static void vectorToCsv(const std::string &_filename, const std::vector<Result> &_data);
    };

    /* Times and overlaps aquired from single test */
    struct SingleTestAcquiredData
    {
        unsigned int overlapped = 0;
        long nanos = 0;
        long overhead = 0;
    };

    /* All times acquired for single strategy from single ball radius */
    struct StrategyAcquiredData {
        CuboidOverlapStrategy * strategy;
        std::vector<double> times;

        StrategyResult generateResult();
    };

    /* All data acquired from single ball radius */
    struct AcquiredData {
        std::string             factoryDesc;
        std::vector<StrategyAcquiredData> strategyDatas;
        std::vector<double>     numOverlapped;
        std::size_t             numAll{};

        AcquiredData(const Context &_context, RSAShapePairFactory &_factory);
        static std::vector<AcquiredData> initVector(const Context &_context, IndependentPairFactory &_factory);
        Result generateResult();
        static std::vector<Result> generateResultVector(std::vector<AcquiredData> &_acquiredDatas);
    };

    void warmUp(RSAShapePairFactory &_factor);
    SingleTestAcquiredData test_single_alg(RSAShapePairFactory &_factory, std::size_t _pairs_to_test);
    void test_single_repeat(const Context &_context, RSAShapePairFactory &_factory, AcquiredData &_acquired_data);

    int main(int argc, char **argv);
}

#endif // _CUBOID_SPEED_TEST_H
