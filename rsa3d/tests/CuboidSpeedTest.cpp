//--------------------------------------------------------------------------------------------
// Comparison of different overlap algorithms speed depending on overlap probability
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include <chrono>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "CuboidSpeedTest.h"


using namespace std::chrono;

namespace
{
    // Struct from single algorithm test
    struct SingleTestResult
    {
        unsigned int overlapped = 0;
        unsigned int nanos = 0;
    };

    // Dummy BoundaryConditions
    class MockBC : public BoundaryConditions
    {
        double distance2(double *p1, double *p2)
        {
            return 0;
        }
        
	    double * getTranslation(double *result, double* p1, double* p2)
	    {
	        result[0] = 0;
	        result[1] = 0;
	        result[2] = 0;
	        return result;
	    }
    };
    
    // Helper method. Performs single test of overlap algorithm (set from the outside)
    //----------------------------------------------------------------------------------------
    SingleTestResult test_single_alg(CuboidPairFactory * _factory, std::size_t _pairs_to_test)
    {
        SingleTestResult result;
        MockBC bc;
        system_clock::time_point time_before, time_after;
        nanoseconds dur;
        
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
	    return result;
    }
}


namespace cube_speedtest
{
    // Returns Quantity computed from vector of samples
    //----------------------------------------------------------------------------------------
    Quantity Quantity::fromSamples(const std::vector<double> & _samples)
    {
        if (_samples.size() == 0)
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
        _stream << _quantity.value << "+-" << _quantity.error;
        return _stream;
    }


    // Performs test for given Cuboid factory and returns results
    //----------------------------------------------------------------------------------------
    CuboidTestData perform(CuboidPairFactory * _factory, std::size_t _pairs_to_test, std::size_t _repeats)
    {
        if (_pairs_to_test == 0)
            throw std::runtime_error("_pairs_to_test == 0");
        if (_repeats == 0)
            throw std::runtime_error("_repeats == 0");
          
        CuboidTestData result;
        SingleTestResult single_result;
        std::vector<double> mine_times, tri_times, SAT_times;
        std::vector<double> num_overlapped;
        
        // Reserve memory for vector of measured times and number of overlaps
        mine_times.reserve(_repeats);
        tri_times.reserve(_repeats);
        SAT_times.reserve(_repeats);
        num_overlapped.reserve(_repeats);
        
        result.numAll = _pairs_to_test;
        
        // Test
        for (std::size_t i = 0; i < _repeats; i++)
        {
            // Test my algorithm
            Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy::MINE);
            single_result = test_single_alg(_factory, _pairs_to_test);
            num_overlapped.push_back(single_result.overlapped);
            mine_times.push_back(single_result.nanos);
            
            // Test triangle algorithm
            Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy::TRI_TRI);
            single_result = test_single_alg(_factory, _pairs_to_test);
            num_overlapped.push_back(single_result.overlapped);
            tri_times.push_back(single_result.nanos);
            
            // Test SAT algorithm
            Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy::SAT);
            single_result = test_single_alg(_factory, _pairs_to_test);
            num_overlapped.push_back(single_result.overlapped);
            SAT_times.push_back(single_result.nanos);
        }
        
        // Calculate and return results
        result.numOverlapped = Quantity::fromSamples(num_overlapped);
        result.mineNs = Quantity::fromSamples(mine_times);
        result.triNs = Quantity::fromSamples(tri_times);
        result.SATNs = Quantity::fromSamples(SAT_times);
         
        return result;
    }


    // Prints results saved in _data onto the standard output
    //----------------------------------------------------------------------------------------
    void print_results(CuboidTestData _data)
    {
         
    }


    // Stores all results in vector to csv file _file
    //----------------------------------------------------------------------------------------
    void to_csv(const std::fstream & _file, const std::vector<CuboidTestData> & _data)
    {

    }
}
