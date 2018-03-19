//--------------------------------------------------------------------------------------------
// Mock BoundaryConditions for testing purpose
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _MOCK_BC_H
    #define _MOCK_BC_H
    
#include "../../rsa3d/BoundaryConditions.h"
    
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

#endif //_MOCK_BC_H
