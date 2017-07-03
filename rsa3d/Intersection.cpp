//--------------------------------------------------------------------------------------------
// Packet of functions for detecting intersections between primitives
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include <cmath>
#include <algorithm>
#include <utility>

#include "Intersection.h"


namespace intersection
{
    // Check intersection between two R3 triangles _tri1 and _tri2
    //--------------------------------------------------------------------------------------------
    bool tri_tri3D(const Vector * _tri1, const Vector * _tri2)
    {
        double d1, d2;
        double dtri1[3], dtri2[3];
        double ptri1[3], ptri2[3];
        double t1_tri1, t2_tri1, t1_tri2, t2_tri2;
        double int_min, int_max;
        
        // _tri2 plane and _tri1 signed distances
        Vector N2 = (_tri2[1] - _tri2[0]) ^ (_tri2[2] - _tri2[0]);
        d2 = -(N2 * _tri2[0]);
        for (int i = 0; i < 3; i++)
            dtri1[i] = N2 * _tri1[i] + d2;
        // Early _tri1 rejection test
        if ((dtri1[0] > 0 && dtri1[1] > 0 && dtri1[2] > 0) || (dtri1[0] < 0 && dtri1[1] < 0 && dtri1[2] < 0))
            return false;
        
        // _tri1 plane and _tri2 signed distances
        Vector N1 = (_tri1[1] - _tri1[0]) ^ (_tri1[2] - _tri1[0]);
        d1 = -(N1 * _tri1[0]);
        for (int i = 0; i < 3; i++)
            dtri2[i] = N1 * _tri2[i] + d1;
        // Early _tri2 rejection test
        if ((dtri2[0] > 0 && dtri2[1] > 0 && dtri2[2] > 0) || (dtri2[0] < 0 && dtri2[1] < 0 && dtri2[2] < 0))
            return false;
        
        // Co-planar line and vertices projection
        Vector D = N1 ^ N2;
        int coord;
        if (std::abs(D(0)) >= std::abs(D(1)) && std::abs(D(0)) >= std::abs(D(2)))
            coord = 0;
        else if (std::abs(D(1)) >= std::abs(D(2)))
            coord = 1;
        else
            coord = 2;
        
        for (int i = 0; i < 3; i++) {
            ptri1[i] = _tri1[i](coord);
            ptri2[i] = _tri2[i](coord);
        }
        
        // Line intervals for both triangles
        t1_tri1 = ptri1[0] + (ptri1[1] - ptri1[0]) * dtri1[0] / (dtri1[0] - dtri1[1]);
        t2_tri1 = ptri1[0] + (ptri1[2] - ptri1[0]) * dtri1[0] / (dtri1[0] - dtri1[2]);
        t1_tri2 = ptri2[0] + (ptri2[1] - ptri2[0]) * dtri2[0] / (dtri2[0] - dtri2[1]);
        t2_tri2 = ptri2[0] + (ptri2[2] - ptri2[0]) * dtri2[0] / (dtri2[0] - dtri2[2]);
        
        // Check interval intersection
        if (t2_tri1 < t1_tri1)
            std::swap(t1_tri1, t2_tri1);
        if (t2_tri2 < t1_tri2)
            std::swap(t1_tri2, t2_tri2);
        int_min = std::max(t1_tri1, t1_tri2);
        int_max = std::min(t2_tri1, t2_tri2);
        if (int_min > int_max)
            return false;
        else
            return true;
    } 
}
