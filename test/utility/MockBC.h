//--------------------------------------------------------------------------------------------
// Mock BoundaryConditions for testing purpose
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _MOCK_BC_H
    #define _MOCK_BC_H
    
#include "../../rsa3d/BoundaryConditions.h"
#include <cstdlib>
    
class MockBC : public BoundaryConditions
{
    double distance2(const double *p1, const double *p2) { return 0; }
    
    double * getTranslation(double *result, const double* p1, const double* p2) {
        for (std::size_t i = 0; i < RSA_SPATIAL_DIMENSION; i++)
            result[i] = 0;
        return result;
    }
};

#endif //_MOCK_BC_H
