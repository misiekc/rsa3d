//--------------------------------------------------------------------------------------------
// Comparison of different overlap algorithms speed depending on overlap probability
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include <chrono>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <bits/unique_ptr.h>
#include <fstream>
#include <iterator>

#include "CuboidSpeedTest.h"
#include "utility/MockBC.h"
#include "../shapes/cube_strategies/MineOverlap.h"
#include "../shapes/cube_strategies/TriTriOverlap.h"
#include "../shapes/cube_strategies/SATOverlap.h"
#include "../shapes/cube_strategies/OptimizedSATOverlap.h"


using namespace std::chrono;

namespace
{
    // Struct for single algorithm test
    struct SingleTestResult
    {
        unsigned int overlapped = 0;
        long nanos = 0;
        long overhead = 0;
    };

    class Timer
    {
        system_clock::time_point startTime;
        system_clock::time_point endTime;

    public:
        void start() {
            startTime = system_clock::now();
        }

        void stop() {
            endTime = system_clock::now();
        }

        template <typename DUR>
        long count() {
            DUR duration = duration_cast<DUR>(endTime - startTime);
            return duration.count();
        }
    };

    // Helper method. Performs single test of overlap algorithm (set from the outside)
    //----------------------------------------------------------------------------------------
    SingleTestResult test_single_alg(CuboidPairFactory * _factory, std::size_t _pairs_to_test)
    {
        SingleTestResult result{};
        OverlapStrategy * strategy = Cuboid::getOverlapStrategy();
        MockBC bc;
        Timer timer;

        std::cout << std::setw(20) << std::left << Cuboid::getOverlapStrategy()->getName() << " test... " << std::flush;

        timer.start();
        for (std::size_t i = 0; i < _pairs_to_test; i++) {
            CuboidPairFactory::CuboidPair pair;
            pair = _factory->generate();
            if (pair.first->overlap(&bc, pair.second))
                result.overlapped++;
            pair.free();
        }
        timer.stop();
        result.nanos = timer.count<nanoseconds>() / _pairs_to_test;

        timer.start();
        for (std::size_t i = 0; i < _pairs_to_test; i++) {
            CuboidPairFactory::CuboidPair pair;
            pair = _factory->generate();
            strategy->runOverheadOperations(pair.first, pair.second);
            pair.free();
        }
        timer.stop();
        result.overhead = timer.count<nanoseconds>() / _pairs_to_test;
        result.nanos -= result.overhead;
	    
	    std::cout << result.nanos << " ns, " << result.overhead << " ns overhead, " << result.overlapped << " overlapped" << std::endl;
	    return result;
    }

    void die(const std::string & reason)
    {
        std::cerr << reason << std::endl;
        exit(EXIT_FAILURE);
    }

    OverlapStrategy * strategyFromString(const std::string &_name) {
        if (_name == "mine")
            return new MineOverlap;
        else if (_name == "sat")
            return new SATOverlap;
        else if (_name == "optimised_sat")
            return new OptimizedSATOverlap;
        else if (_name == "tri_tri")
            return new TriTriOverlap;
        else
            throw std::runtime_error("unknown strategy: " + _name);
    }
}


namespace cube_speedtest
{
    // Perforsm one warm up test for consistent results. The first one seems to give slower
    // result regardless of chosen algorithm.
    //----------------------------------------------------------------------------------------
    void warmUp(CuboidPairFactory * _factory)
    {
        std::cout << "Warm up test..." << std::endl;
        test_single_alg(_factory, 200000);
    }


    // Performs test for given Cuboid factory and returns results
    //----------------------------------------------------------------------------------------
    void test_single_repeat(CuboidPairFactory *_factory, AcquiredData &_acquired_data,
                            const std::vector<OverlapStrategy *> &_strategies, std::size_t _pairs_to_test)
    {
        if (_pairs_to_test == 0)
            throw std::runtime_error("_pairs_to_test == 0");

        std::cout << std::endl;
//        std::cout << "Starting (" << _factory->getDescription() << ", pairs: " << _pairs_to_test
//            << ", repeats: " << _repeats << ")..." << std::endl;

        /*std::stringstream tryNoInfoStream;
        tryNoInfoStream << "(" << (i + 1) << "/" << _repeats << ") ";
        std::string tryNoInfo = tryNoInfoStream.str();
        std::string tryNoInfoSpace(tryNoInfo.length(), ' ');*/

        // Test each strategy
        for (std::size_t j = 0; j < _strategies.size(); j++) {
            Cuboid::setOverlapStrategy(_strategies[j]);
//            std::cout << (j == 0 ? tryNoInfo : tryNoInfoSpace);
            SingleTestResult single_result = test_single_alg(_factory, _pairs_to_test);

            _acquired_data.strategyDatas[j].times.push_back(single_result.nanos);
            _acquired_data.numOverlapped.push_back(single_result.overlapped);
        }
    }


    // Prints results saved in _data onto the standard output
    //----------------------------------------------------------------------------------------
    void Result::printResults()
    {
        std::cout << "(" << this->numOverlapped << ")/" << this->numAll << " overlapped. Overlap probability: "
            << this->overlapProb << std::endl;
        for (auto strategyData : strategyResults)
            std::cout << std::left << std::setw(30) << (strategyData.strategy->getName() + " avg. time") << ": " << strategyData.time << std::endl;
        std::cout << std::left << std::setw(30) << "Factory used" << ": " << this->factoryDesc << std::endl;
    }


    // Stores all results in vector to csv output stream _out
    //----------------------------------------------------------------------------------------
    void Result::toCsv(std::ostream & _out, const std::vector<Result> & _data)
    {
        _out << "probability,error";
        for (auto strategyData : _data[0].strategyResults)
            _out << "," << strategyData.strategy->getName() << " time," << strategyData.strategy->getName() << " error";
        _out << ",num of overlapped,error,num of all,factory desc" << std::endl;
            
        for (auto d : _data) {
            _out << d.overlapProb.value << ",";
            _out << d.overlapProb.error << ",";
            for (auto strategyData : d.strategyResults) {
                _out << strategyData.time.value << ",";
                _out << strategyData.time.error << ",";
            }
            _out << d.numOverlapped.value << ",";
            _out << d.numOverlapped.error << ",";
            _out << d.numAll << ",";
            _out << d.factoryDesc << std::endl;
        }
        
        _out.flush();
    }

    // Performs cuboid speedtests with parameters passed to process
    //----------------------------------------------------------------------------------------
    int main(int argc, char **argv) {
        if (argc < 3)
            die("Usage: ./rsa cube_speedtest [input] [output = cube_speedtest.csv]");

        std::ifstream input(argv[2]);
        if (!input)
            die("Error opening " + std::string(argv[2]) + " file to read");

        auto config = std::unique_ptr<Config> (Config::parse(input));
        input.close();

        std::string output = "cube_speedtest.csv";
        if (argc == 4)
            output = argv[3];

        std::size_t pairs = config->getUnsignedInt("pairs");
        std::size_t repeats = config->getUnsignedInt("repeats");

        ShapeFactory::initShapeClass("Cuboid", "3 " + config->getString("cuboid_size"));
        BallFactory * factory = BallFactory::getInstance();

        std::istringstream strategiesStream(config->getString("strategies"));
        std::string strategyName;
        std::vector<OverlapStrategy *> strategies;
        while (strategiesStream >> strategyName)
            strategies.push_back(strategyFromString(strategyName));

        std::istringstream ballRadiaStream(config->getString("ball_radia"));
        std::vector<double> ballRadia{std::istream_iterator<double>(ballRadiaStream), std::istream_iterator<double>()};

        std::vector<AcquiredData> acquiredDatas;
        std::transform(ballRadia.begin(), ballRadia.end(), std::back_inserter(acquiredDatas),
                       [&](double radius) {
                           factory->setRadius(radius);
                           return AcquiredData(factory, strategies, pairs);
                       });

        // Warm up and perform tests
        warmUp(factory);
        for (size_t i = 0; i < repeats; i++) {
            for (size_t j = 0; j < acquiredDatas.size(); j++) {
                factory->setRadius(ballRadia[j]);
                test_single_repeat(factory, acquiredDatas[j], strategies, pairs);
            }
        }

        std::vector<Result> results{};
        std::transform(acquiredDatas.begin(), acquiredDatas.end(), std::back_inserter(results),
                       [](AcquiredData acquiredData) {
                           return acquiredData.generateResult();
                       });

        // Print results
        std::cout << std::endl;
        std::cout << ">> Obtained results:" << std::endl << std::endl;
        for (auto data : results) {
            data.printResults();
            std::cout << std::endl;
        }

        // Print aquired probabilities
        std::cout << ">> Probabilities for ball radia" << std::endl;
        std::cout << std::left;
        for (auto data : results) {
            std::cout << std::setw(15) << data.overlapProb.value << " : " << data.factoryDesc << std::endl;
        }
        std::cout << std::endl;

        // Store to file
        std::cout << ">> Storing to file " << output << "..." << std::endl;
        std::ofstream file(output);
        if (!file)
            die("Error opening " + output + " file to write");

        cube_speedtest::Result::toCsv(file, results);
        file.close();

        for (auto strategy : strategies)
            delete strategy;

        return EXIT_SUCCESS;
    }

    StrategyResult StrategyAcquiredData::generateResult() {
        return StrategyResult{strategy, Quantity::fromSamples(times)};
    }

    Result AcquiredData::generateResult() {
        Result result;
        result.numAll = this->numAll;
        result.factoryDesc = this->factoryDesc;
        result.numOverlapped = Quantity::fromSamples(this->numOverlapped);
        result.overlapProb = Quantity(result.numOverlapped.value / this->numAll,
                                      result.numOverlapped.error / this->numAll);

        for (auto strategyAcquiredData : this->strategyDatas)
            result.strategyResults.push_back(strategyAcquiredData.generateResult());

        return result;
    }
}
