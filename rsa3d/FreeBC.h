//--------------------------------------------------------------------------------------------
// Mock BoundaryConditions for testing purpose
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _MOCK_BC_H
    #define _MOCK_BC_H
    
#include "BoundaryConditions.h"
#include <cstdlib>

template <unsigned short DIMENSION>
class FreeBC : public BoundaryConditions<DIMENSION>
{
public:
    using vector = Vector<DIMENSION>;

    double distance2(const vector &p1, const vector &p2) const { return (p1 - p2).norm2(); }
    vector getTranslation(const vector &p1, const vector &p2) const { return vector(); }
};

using RSAFreeBC = FreeBC<RSA_SPATIAL_DIMENSION>;

#endif //_MOCK_BC_H
