//--------------------------------------------------------------------------------------------
// Comparison of different overlap algorithms speed depending on overlap probability
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include "CuboidPairFactory.h"

#include <vector>
#include <fstream>


namespace cube_speedtest 
{   
    // Struct representing quantity with uncertainity
    template
    <typename NUM>
    struct Quantity {
        NUM value;
        NUM error;
        
        Quantity fromSamples(const std::vector<NUM> & _samples);
    };

    // Struct representing analysis results
    struct CuboidTestData {
        std::size_t     numOverlapped;
        std::size_t     numAll;
        float           overlapProb;
        Quantity<float> mineNs;
        Quantity<float> triNs;
        Quantity<float> SATNs;
    };

    CuboidTestData perform(CuboidPairFactory * _factory, std::size_t _pairs_to_test, std::size_t _repeats);
    void print_results(CuboidTestData _data);
    void to_csv(const std::fstream & _file, const std::vector<CuboidTestData> & _data);   
}
