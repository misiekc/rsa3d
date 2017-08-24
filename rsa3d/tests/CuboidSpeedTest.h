//--------------------------------------------------------------------------------------------
// Comparison of different overlap algorithms speed depending on overlap probability
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include "CuboidPairFactory.h"

#include <vector>
#include <fstream>
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
    struct CuboidTestData {
        Quantity        numOverlapped;
        std::size_t     numAll;
        Quantity        mineNs;
        Quantity        triNs;
        Quantity        SATNs;
    };

    CuboidTestData perform(CuboidPairFactory * _factory, std::size_t _pairs_to_test, std::size_t _repeats);
    void print_results(CuboidTestData _data);
    void to_csv(const std::fstream & _file, const std::vector<CuboidTestData> & _data);   
}
