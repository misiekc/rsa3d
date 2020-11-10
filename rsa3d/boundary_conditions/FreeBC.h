//--------------------------------------------------------------------------------------------
// Mock BoundaryConditions for testing purpose
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef BOUNDARY_CONDITIONS_FREEBC_H_
    #define BOUNDARY_CONDITIONS_FREEBC_H_
    
#include "../BoundaryConditions.h"
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

#endif //BOUNDARY_CONDITIONS_FREEBC_H_
