//--------------------------------------------------------------------------------------------
// Comparison of different overlap algorithms speed depending on overlap probability
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CUBOID_SPEED_TEST_H
    #define _CUBOID_SPEED_TEST_H


#include "../Config.h"
#include "../ShapeFactory.h"
#include "utility/BallFactory.h"
#include "utility/CuboidPairFactory.h"
#include "utility/Quantity.h"

#include <vector>
#include <ostream>


namespace cube_speedtest
{
    /* Test context */
    struct Context {
        std::size_t pairs;
        std::size_t repeats;
        std::vector<OverlapStrategy *> strategies;
        std::vector<double> ballRadia;
        BallFactory * factory;

        void load(std::istream & _input);
    };

    /* End results for single strategy */
    struct StrategyResult {
        OverlapStrategy * strategy;
        Quantity        time;
    };

    /* End results for single pair factory size */
    struct Result {
        Quantity        numOverlapped;
        std::size_t     numAll;
        Quantity        overlapProb;
        std::vector<StrategyResult>   strategyResults;
        std::string     factoryDesc;
        
        void printResults();
        static void toCsv(std::ostream & _out, const std::vector<Result> & _data);
    };

    /* Times acquired from single strategy */
    struct StrategyAcquiredData {
        OverlapStrategy * strategy;
        std::vector<double> times;

        StrategyResult generateResult();
    };

    /* Data acquired from single pair factory size */
    struct AcquiredData {
        std::string             factoryDesc;
        std::vector<StrategyAcquiredData> strategyDatas;
        std::vector<double>     numOverlapped;
        std::size_t             numAll;

        AcquiredData(CuboidPairFactory * _factory, std::vector<OverlapStrategy *> _strategies, std::size_t _pairs_to_test) {
            this->factoryDesc = _factory->getDescription();
            this->numAll = _pairs_to_test;
            for (auto strategy : _strategies)
                this->strategyDatas.push_back(StrategyAcquiredData{strategy});
        }
        Result generateResult();
    };

    void warmUp(CuboidPairFactory * _factor);
    void test_single_repeat(CuboidPairFactory *_factory, AcquiredData &_acquired_data,
                                const std::vector<OverlapStrategy *> &_strategies, std::size_t _pairs_to_test);

    int main(int argc, char **argv);
}

#endif // _CUBOID_SPEED_TEST_H
