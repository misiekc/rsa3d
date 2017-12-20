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
    };
    
    // Helper method. Performs single test of overlap algorithm (set from the outside)
    //----------------------------------------------------------------------------------------
    SingleTestResult test_single_alg(CuboidPairFactory * _factory, std::size_t _pairs_to_test)
    {
        SingleTestResult result{};
        MockBC bc;
        system_clock::time_point time_before, time_after;
        nanoseconds dur;

        std::cout << std::setw(20) << std::left << Cuboid::getOverlapStrategy()->getName()
                  << " test... " << std::flush;

        time_before = system_clock::now();
	    for (std::size_t i = 0; i < _pairs_to_test; i++) {
	        CuboidPairFactory::CuboidPair pair;
	        pair = _factory->generate();
	        if (pair.first->overlap(&bc, pair.second))
	            result.overlapped++;
	        pair.free();
	    }
	    time_after = system_clock::now();
	    dur = duration_cast<std::chrono::nanoseconds>(time_after - time_before);
	    dur /= _pairs_to_test;
	    result.nanos = dur.count();
	    
	    std::cout << result.nanos << " ns, " << result.overlapped << " overlapped" << std::endl;
	    return result;
    }
}


namespace cube_speedtest
{
    // Returns Quantity computed from vector of samples
    //----------------------------------------------------------------------------------------
    Quantity Quantity::fromSamples(const std::vector<double> & _samples)
    {
        if (_samples.empty())
            return Quantity();
        else if (_samples.size() == 1)
            return Quantity(_samples[0], 0);
            
        Quantity result;
        double sum, dev_sum;
        
        sum = std::accumulate(_samples.begin(), _samples.end(), 0);
        result.value = sum / _samples.size();
        dev_sum = std::accumulate(_samples.begin(), _samples.end(), 0,
            [&](double s, double next) {
                s += pow(result.value - next, 2);
                return s;
            });
       result.error = dev_sum / _samples.size() / (_samples.size() - 1);
       result.error = sqrt(result.error);
       return result;
    }
    

    // Friend stream insertion operator
    //----------------------------------------------------------------------------------------
    std::ostream & operator<<(std::ostream & _stream, const Quantity & _quantity)
    {
        _stream << _quantity.value << " +- " << _quantity.error;
        return _stream;
    }

    
    // Perforsm one warm up test for consistent results. The first one seems to give slower
    // result regardless of chosen algorithm.
    //----------------------------------------------------------------------------------------
    void warmUp(CuboidPairFactory * _factory)
    {
        std::cout << "Warm up test..." << std::endl;
        test_single_alg(_factory, 2000000);
    }


    // Performs test for given Cuboid factory and returns results
    //----------------------------------------------------------------------------------------
    TestData perform(CuboidPairFactory *_factory, const std::vector<OverlapStrategy *> &_strategies, std::size_t _pairs_to_test,
                     std::size_t _repeats)
    {
        if (_pairs_to_test == 0)
            throw std::runtime_error("_pairs_to_test == 0");
        if (_repeats == 0)
            throw std::runtime_error("_repeats == 0");

        std::cout << std::endl;
        std::cout << "Starting (" << _factory->getDescription() << ", pairs: " << _pairs_to_test
            << ", repeats: " << _repeats << ")..." << std::endl;

        // Test
        std::vector<double> num_overlapped;
        std::vector<std::vector<double>> times(_strategies.size());
        for (std::size_t i = 0; i < _repeats; i++)
        {
            std::stringstream tryNoInfoStream;
            tryNoInfoStream << "(" << (i + 1) << "/" << _repeats << ") ";
            std::string tryNoInfo = tryNoInfoStream.str();
            std::string tryNoInfoSpace(tryNoInfo.length(), ' ');

            // Test each strategy
            for (std::size_t j = 0; j < _strategies.size(); j++) {
                Cuboid::setOverlapStrategy(_strategies[j]);
                std::cout << (j == 0 ? tryNoInfo : tryNoInfoSpace);
                SingleTestResult single_result = test_single_alg(_factory, _pairs_to_test);
                num_overlapped.push_back(single_result.overlapped);
                times[j].push_back(single_result.nanos);
            }
        }

        // Calculate and return results
        TestData result;
        result.numAll = _pairs_to_test;
        result.factoryDesc = _factory->getDescription();
        result.numOverlapped = Quantity::fromSamples(num_overlapped);
        result.overlapProb = Quantity(result.numOverlapped.value / _pairs_to_test,
                                      result.numOverlapped.error / _pairs_to_test);

        for (std::size_t i = 0; i < _strategies.size(); i++) {
            result.strategyDatas.push_back(StrategyData {
                    _strategies[i]->getName(),
                    Quantity::fromSamples(times[i]) });
        }

        return result;
    }


    // Prints results saved in _data onto the standard output
    //----------------------------------------------------------------------------------------
    void TestData::printResults()
    {
        std::cout << "(" << this->numOverlapped << ")/" << this->numAll << " overlapped. Overlap probability: "
            << this->overlapProb << std::endl;
        for (auto strategyData : strategyDatas)
            std::cout << std::left << std::setw(30) << (strategyData.name + " avg. time") << ": " << strategyData.time << std::endl;
        std::cout << std::left << std::setw(30) << "Factory used" << ": " << this->factoryDesc << std::endl;
    }


    // Stores all results in vector to csv output stream _out
    //----------------------------------------------------------------------------------------
    void TestData::toCsv(std::ostream & _out, const std::vector<TestData> & _data)
    {
        _out << "probability,error";
        for (auto strategyData : _data[0].strategyDatas)
            _out << "," << strategyData.name << " time," << strategyData.name << " error";
        _out << ",num of overlapped,error,num of all,factory desc" << std::endl;
            
        for (auto d : _data) {
            _out << d.overlapProb.value << ",";
            _out << d.overlapProb.error << ",";
            for (auto strategyData : d.strategyDatas) {
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
}
