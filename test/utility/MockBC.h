//--------------------------------------------------------------------------------------------
// Mock BoundaryConditions for testing purpose
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _MOCK_BC_H
    #define _MOCK_BC_H
    
#include "../../rsa3d/BoundaryConditions.h"
#include <cstdlib>

template <unsigned short DIMENSION>
class MockBC : public BoundaryConditions<DIMENSION>
{
public:
    using vector = Vector<DIMENSION>;

    double distance2(const vector &p1, const vector &p2) const { return 0; }
    vector getTranslation(const vector &p1, const vector &p2) const { return vector(); }
};

using RSAMockBC = MockBC<RSA_SPATIAL_DIMENSION>;

#endif //_MOCK_BC_H
