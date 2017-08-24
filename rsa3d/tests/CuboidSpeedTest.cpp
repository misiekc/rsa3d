//--------------------------------------------------------------------------------------------
// Comparison of different overlap algorithms speed depending on overlap probability
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include "CuboidSpeedTest.h"


namespace
{
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
}


namespace cube_speedtest
{
    // Returns Quantity computed from vector of samples
    //----------------------------------------------------------------------------------------
    template
    <typename NUM>
    Quantity<NUM> Quantity<NUM>::fromSamples(const std::vector<NUM> & _samples)
    {
        return cube_speedtest::Quantity<NUM>();
    }


    // Perform test for given Cuboid factory and return results
    //----------------------------------------------------------------------------------------
    CuboidTestData perform(CuboidPairFactory * _factory, std::size_t _pairs_to_test, std::size_t _repeats)
    {
        return CuboidTestData();
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
