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
#include "MockBC.h"
#include "../shapes/cube_strategies/MineOverlap.h"


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

        std::cout << std::setw(15) << std::left << Cuboid::getOverlapStrategy()->getName()
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
    TestData perform(CuboidPairFactory * _factory, std::size_t _pairs_to_test, std::size_t _repeats)
    {
        if (_pairs_to_test == 0)
            throw std::runtime_error("_pairs_to_test == 0");
        if (_repeats == 0)
            throw std::runtime_error("_repeats == 0");
          
        TestData result;
        SingleTestResult single_result;
        std::vector<double> mine_times, tri_times, SAT_times;
        std::vector<double> num_overlapped;
        
        // Reserve memory for vector of measured times and number of overlaps
        mine_times.reserve(_repeats);
        tri_times.reserve(_repeats);
        SAT_times.reserve(_repeats);
        num_overlapped.reserve(_repeats);
        
        result.numAll = _pairs_to_test;
        result.factoryDesc = _factory->getDescription();
        std::cout << std::endl;
        std::cout << "Starting (" << result.factoryDesc << ", pairs: " << _pairs_to_test
            << ", repeats: " << _repeats << ")..." << std::endl;

        MineOverlap mineOverlap;

        // Test
        for (std::size_t i = 0; i < _repeats; i++)
        {
            std::stringstream tryNoInfoStream;
            tryNoInfoStream << "(" << (i + 1) << "/" << _repeats << ") ";
            std::string tryNoInfo = tryNoInfoStream.str();
            std::string tryNoInfoSpace(tryNoInfo.length(), ' ');
            
            // Test my algorithm
            Cuboid::setOverlapStrategy(&mineOverlap);
            std::cout << tryNoInfo;
            single_result = test_single_alg(_factory, _pairs_to_test);
            num_overlapped.push_back(single_result.overlapped);
            mine_times.push_back(single_result.nanos);
            
            // Test triangle algorithm
            Cuboid::setOverlapStrategy(&mineOverlap);
            std::cout << tryNoInfoSpace;
            single_result = test_single_alg(_factory, _pairs_to_test);
            num_overlapped.push_back(single_result.overlapped);
            tri_times.push_back(single_result.nanos);
            
            // Test SAT algorithm
            Cuboid::setOverlapStrategy(&mineOverlap);
            std::cout << tryNoInfoSpace;
            single_result = test_single_alg(_factory, _pairs_to_test);
            num_overlapped.push_back(single_result.overlapped);
            SAT_times.push_back(single_result.nanos);
        }
        
        // Calculate and return results
        result.numOverlapped = Quantity::fromSamples(num_overlapped);
        result.mineNs = Quantity::fromSamples(mine_times);
        result.triNs = Quantity::fromSamples(tri_times);
        result.SATNs = Quantity::fromSamples(SAT_times);
        result.overlapProb = Quantity(result.numOverlapped.value / _pairs_to_test, 
            result.numOverlapped.error / _pairs_to_test);
         
        return result;
    }


    // Prints results saved in _data onto the standard output
    //----------------------------------------------------------------------------------------
    void TestData::printResults()
    {
        std::cout << "(" << this->numOverlapped << ")/" << this->numAll << " overlapped. Overlap probability: "
            << this->overlapProb << std::endl;
        std::cout << "Mine avg. time    : " << this->mineNs << std::endl;
        std::cout << "Tri-tri avg. time : " << this->triNs << std::endl;
        std::cout << "SAT avg. time     : " << this->SATNs << std::endl;
        std::cout << "Factory used      : " << this->factoryDesc << std::endl;
    }


    // Stores all results in vector to csv output stream _out
    //----------------------------------------------------------------------------------------
    void TestData::toCsv(std::ostream & _out, const std::vector<TestData> & _data)
    {
        _out << "probability,error,my alg time,error,tri-tri alg time,error,SAT alg time,error,"
            "num of overlapped,error,num of all,factory desc" << std::endl;
            
        for (auto d : _data) {
            _out << d.overlapProb.value << ",";
            _out << d.overlapProb.error << ",";
            _out << d.mineNs.value << ",";
            _out << d.mineNs.error << ",";
            _out << d.triNs.value << ",";
            _out << d.triNs.error << ",";
            _out << d.SATNs.value << ",";
            _out << d.SATNs.error << ",";
            _out << d.numOverlapped.value << ",";
            _out << d.numOverlapped.error << ",";
            _out << d.numAll << ",";
            _out << d.factoryDesc << std::endl;
        }
        
        _out.flush();
    }
}
