//--------------------------------------------------------------------------------------------
// Comparison of different overlap algorithms speed depending on overlap probability
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <bits/unique_ptr.h>
#include <fstream>
#include <iterator>

#include "CuboidSpeedTest.h"
#include "utility/MockBC.h"
#include "../rsa3d/shape/shapes/cuboid/MineOverlap.h"
#include "../rsa3d/shape/shapes/cuboid/TriTriOverlap.h"
#include "../rsa3d/shape/shapes/cuboid/SATOverlap.h"
#include "../rsa3d/shape/shapes/cuboid/OptimizedSATOverlap.h"
#include "../rsa3d/Timer.h"
#include "../rsa3d/Utils.h"

using namespace std::chrono;

namespace
{

    CuboidOverlapStrategy * strategyFromString(const std::string &_name) {
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
    // Performs one warm up test for consistent results. The first one seems to claim slower
    // speed regardless of chosen algorithm.
    //----------------------------------------------------------------------------------------
    void warmUp(RSAShapePairFactory &_factory)
    {
        std::cout << "Warm up test..." << std::endl;
        test_single_alg(_factory, 2000000);
    }

    // Performs single test of overlap algorithm (set from the outside)
    //----------------------------------------------------------------------------------------
    SingleTestAcquiredData test_single_alg(RSAShapePairFactory &_factory, std::size_t _pairs_to_test)
    {
        SingleTestAcquiredData result{};
        CuboidOverlapStrategy * strategy = Cuboid::getOverlapStrategy();
        RSAMockBC bc;
        Timer timer;

        std::cout << std::setw(20) << std::left << Cuboid::getOverlapStrategy()->getName() << " test... " << std::flush;

        timer.start();
        for (std::size_t i = 0; i < _pairs_to_test; i++) {
            RSAShapePairFactory::ShapePair pair = _factory.generate();
            if (pair.first()->overlap(&bc, pair.second()))
                result.overlapped++;
        }
        timer.stop();
        result.nanos = timer.count<nanoseconds>() / _pairs_to_test;

        timer.start();
        for (std::size_t i = 0; i < _pairs_to_test; i++) {
            RSAShapePairFactory::ShapePair pair = _factory.generate();
            strategy->runOverheadOperations((Cuboid *)pair.first(), (Cuboid *)pair.second());
        }
        timer.stop();
        result.overhead = timer.count<nanoseconds>() / _pairs_to_test;
        result.nanos -= result.overhead;

        std::cout << result.nanos << " ns, " << result.overhead << " ns overhead, " << result.overlapped << " overlapped" << std::endl;
        return result;
    }


    // Performs single repeat of _pairs_to_test test using given _factory for all passed
    // _strategies and stored result in _acquired_data
    //----------------------------------------------------------------------------------------
    void test_single_repeat(const Context &_context, RSAShapePairFactory &_factory, AcquiredData &_acquired_data)
    {
        // Test each strategy
        for (std::size_t j = 0; j < _context.strategies.size(); j++) {
            Cuboid::setOverlapStrategy(_context.strategies[j]);
            SingleTestAcquiredData single_result = test_single_alg(_factory, _context.pairs);

            _acquired_data.strategyDatas[j].times.push_back(single_result.nanos);
            _acquired_data.numOverlapped.push_back(single_result.overlapped);
        }
    }


    // Performs cuboid speedtests with parameters passed to process
    //----------------------------------------------------------------------------------------
    int main(int argc, char **argv) {
        if (argc < 3)
            die("Usage: ./rsa cube_speedtest [input] [output = cube_speedtest.csv]");

        Context context;
        context.load(argv[2]);
        BallFactory factory;
        std::vector<AcquiredData> acquiredDatas = AcquiredData::initVector(context, factory);

        // Warm up and perform tests
        warmUp(factory);
        std::cout << std::endl;
        for (size_t i = 0; i < context.repeats; i++) {
            for (size_t j = 0; j < context.ballRadia.size(); j++) {
                factory.setRadius(context.ballRadia[j]);
                std::cout << std::endl << ">> Repeat " << (i + 1) << "/" << context.repeats
                          << " for " << factory.getDescription() << "..." << std::endl;
                test_single_repeat(context, factory, acquiredDatas[j]);
            }
            std::cout << std::endl;
        }

        // Print results on std out
        std::vector<Result> results = AcquiredData::generateResultVector(acquiredDatas);
        Result::printVector(results);

        // Store to file
        if (argc == 4)
            cube_speedtest::Result::vectorToCsv(argv[3], results);
        else
            cube_speedtest::Result::vectorToCsv("cube_speedtest.csv", results);
        return EXIT_SUCCESS;
    }


    // Load test context from file _filename
    //----------------------------------------------------------------------------------------
    void Context::load(const std::string &_filename) {
        std::ifstream input(_filename);
        if (!input)
            die("Error opening " + std::string(_filename) + " file to read");
        auto config = Config::parse(input);

        this->pairs = config->getUnsignedInt("pairs");
        this->repeats = config->getUnsignedInt("repeats");

        // Load strategies
        std::istringstream strategiesStream(config->getString("strategies"));
        std::string strategyName;
        while (strategiesStream >> strategyName)
            this->strategies.push_back(strategyFromString(strategyName));

        // Load ball radia
        std::istringstream ballRadiaStream(config->getString("ball_radia"));
        std::copy(std::istream_iterator<double>(ballRadiaStream),
                  std::istream_iterator<double>(),
                  std::back_inserter(this->ballRadia));

        ShapeFactory::initShapeClass("Cuboid", "3 " + config->getString("cuboid_size"));
        delete config;
    }

    Context::~Context() {
        for (auto strategy : this->strategies)
            delete strategy;
    }


    // Prints results saved in _data onto the standard output
    //----------------------------------------------------------------------------------------
    void Result::print() const
    {
        std::cout << "(" << this->numOverlapped << ")/" << this->numAll << " overlapped. Overlap probability: "
            << this->overlapProb << std::endl;
        for (auto strategyData : strategyResults)
            std::cout << std::left << std::setw(30) << (strategyData.strategy->getName() + " avg. time") << ": " << strategyData.time << std::endl;
        std::cout << std::left << std::setw(30) << "Factory used" << ": " << this->factoryDesc << std::endl;
    }

    // Prints whole vector of Result-s onto standard output
    //----------------------------------------------------------------------------------------
    void Result::printVector(const std::vector<Result> &_results) {
        std::cout << std::endl;
        std::cout << ">> Obtained results:" << std::endl << std::endl;
        for (const auto &data : _results) {
            data.print();
            std::cout << std::endl;
        }

        // Print aquired probabilities
        std::cout << ">> Probabilities for ball radia" << std::endl;
        std::cout << std::left;
        for (auto data : _results)
            std::cout << std::setw(15) << data.overlapProb.value << " : " << data.factoryDesc << std::endl;
        std::cout << std::endl;
    }

    // Stores all results in vector to csv output stream _out
    //----------------------------------------------------------------------------------------
    void Result::vectorToCsv(const std::string &_filename, const std::vector<Result> &_data)
    {
        std::ofstream file(_filename);
        if (!file)
            die("Error opening " + _filename + " file to write");
        std::cout << ">> Storing to file " << _filename << "..." << std::endl;

        file << "probability,error";
        for (auto strategyData : _data[0].strategyResults)
            file << "," << strategyData.strategy->getName() << " time," << strategyData.strategy->getName() << " error";
        file << ",num of overlapped,error,num of all,factory desc" << std::endl;

        for (auto d : _data) {
            file << d.overlapProb.value << ",";
            file << d.overlapProb.error << ",";
            for (auto strategyData : d.strategyResults) {
                file << strategyData.time.value << ",";
                file << strategyData.time.error << ",";
            }
            file << d.numOverlapped.value << ",";
            file << d.numOverlapped.error << ",";
            file << d.numAll << ",";
            file << d.factoryDesc << std::endl;
        }

        file.close();
    }

    // Inits and returnes a vector for storing data based on _context
    //----------------------------------------------------------------------------------------
    std::vector<AcquiredData> AcquiredData::initVector(const Context &_context, BallFactory &_factory) {
        std::vector<AcquiredData> acquiredDatas{};
        std::transform(_context.ballRadia.begin(), _context.ballRadia.end(), std::back_inserter(acquiredDatas),
                  [&](double radius) {
                      _factory.setRadius(radius);
                      return AcquiredData(_context, _factory);
                  });
        return acquiredDatas;
    }

    // Constructs itself based on given _context
    //----------------------------------------------------------------------------------------
    AcquiredData::AcquiredData(const Context &_context, RSAShapePairFactory &_factory) {
        this->factoryDesc = _factory.getDescription();
        this->numAll = _context.pairs;
        for (auto strategy : _context.strategies)
            this->strategyDatas.push_back(StrategyAcquiredData{strategy});
    }

    // Generates Result from collected data - calculates means and std devs of times and
    // number of overlaps
    //----------------------------------------------------------------------------------------
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

    // Generate vector of Result-s from all passed _acquiredDatas
    //----------------------------------------------------------------------------------------
    std::vector<Result> AcquiredData::generateResultVector(std::vector<AcquiredData> &_acquiredDatas) {
        std::vector<Result> results{};
        std::transform(_acquiredDatas.begin(), _acquiredDatas.end(), std::back_inserter(results),
                       [](AcquiredData acquiredData) {
                           return acquiredData.generateResult();
                       });
        return results;
    }

    // Generates StrategyResult - it calculates mean and std dev of times
    //----------------------------------------------------------------------------------------
    StrategyResult StrategyAcquiredData::generateResult() {
        return StrategyResult{strategy, Quantity::fromSamples(times)};
    }
}
