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

        Result generateResult();
    };

    void warmUp(CuboidPairFactory * _factor);
    Result perform(CuboidPairFactory *_factory, const std::vector<OverlapStrategy *> & _strategies, std::size_t _pairs_to_test,
                         std::size_t _repeats);

    int main(int argc, char **argv);
}

#endif // _CUBOID_SPEED_TEST_H
