//--------------------------------------------------------------------------------------------
// Comparison of different overlap algorithms speed depending on overlap probability
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CUBOID_SPEED_TEST_H
    #define _CUBOID_SPEED_TEST_H
    

#include "CuboidPairFactory.h"

#include <vector>
#include <ostream>


namespace cube_speedtest 
{   
    // Struct representing quantity with uncertainity
    struct Quantity {
        double value = 0;
        double error = 0;
        
        Quantity()
        { }
        
        Quantity(double _value, double _error) : value(_value), error(_error)
        { }
        
        static Quantity fromSamples(const std::vector<double> & _samples);
        friend std::ostream & operator<<(std::ostream & _stream, const Quantity & _quantity);
    };

    // Struct representing analysis results
    struct TestData {
        Quantity        numOverlapped;
        std::size_t     numAll;
        Quantity        overlapProb;
        Quantity        mineNs;
        Quantity        triNs;
        Quantity        SATNs;
        std::string     factoryDesc;
        
        void printResults();
        static void toCsv(std::ostream & _out, const std::vector<TestData> & _data);
    };

    void warmUp(CuboidPairFactory * _factor);
    TestData perform(CuboidPairFactory * _factory, std::size_t _pairs_to_test, std::size_t _repeats);
}

#endif // _CUBOID_SPEED_TEST_H
