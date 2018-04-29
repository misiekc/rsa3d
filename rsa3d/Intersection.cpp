//--------------------------------------------------------------------------------------------
// Packet of functions for detecting intersections between primitives
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include <cmath>
#include <algorithm>
#include <utility>
#include <cassert>
#include <iostream>

#include "Intersection.h"


namespace intersection
{
    // Check intersection between two R3 triangles _tri1 and _tri2
    //--------------------------------------------------------------------------------------------
    bool tri_tri3D(const Vector<3> * _tri1, const Vector<3> * _tri2)
    {
        double d1, d2;
        double dtri1[3], dtri2[3];
        double ptri1[3], ptri2[3];
        double t1_tri1, t2_tri1, t1_tri2, t2_tri2;
        double int_min, int_max;
        
        // _tri2 plane and _tri1 signed distances
        Vector<3> N2 = (_tri2[1] - _tri2[0]) ^ (_tri2[2] - _tri2[0]);
        d2 = -(N2 * _tri2[0]);
        for (int i = 0; i < 3; i++)
            dtri1[i] = N2 * _tri1[i] + d2;
        // Early _tri1 rejection test
        if ((dtri1[0] > 0 && dtri1[1] > 0 && dtri1[2] > 0) || (dtri1[0] < 0 && dtri1[1] < 0 && dtri1[2] < 0))
            return false;
        
        // _tri1 plane and _tri2 signed distances
        Vector<3> N1 = (_tri1[1] - _tri1[0]) ^ (_tri1[2] - _tri1[0]);
        d1 = -(N1 * _tri1[0]);
        for (int i = 0; i < 3; i++)
            dtri2[i] = N1 * _tri2[i] + d1;
        // Early _tri2 rejection test
        if ((dtri2[0] > 0 && dtri2[1] > 0 && dtri2[2] > 0) || (dtri2[0] < 0 && dtri2[1] < 0 && dtri2[2] < 0))
            return false;
        
        // Co-planar line and vertices projection
        Vector<3> D = N1 ^ N2;
        int coord;
        if (std::abs(D[0]) >= std::abs(D[1]) && std::abs(D[0]) >= std::abs(D[2]))
            coord = 0;
        else if (std::abs(D[1]) >= std::abs(D[2]))
            coord = 1;
        else
            coord = 2;
        
        for (int i = 0; i < 3; i++) {
            ptri1[i] = _tri1[i][coord];
            ptri2[i] = _tri2[i][coord];
        }
        
        // Assure that 0-indexed variables represent vertex lying of the opposi side of D
        if ((dtri1[0] > 0 && dtri1[1] > 0) || (dtri1[0] < 0 && dtri1[1] < 0)) {
            std::swap(dtri1[0], dtri1[2]);
            std::swap(ptri1[0], ptri1[2]);
        } else if ((dtri1[0] > 0 && dtri1[2] > 0) || (dtri1[0] < 0 && dtri1[2] < 0)) {
            std::swap(dtri1[0], dtri1[1]);
            std::swap(ptri1[0], ptri1[1]);
        }
        if ((dtri2[0] > 0 && dtri2[1] > 0) || (dtri2[0] < 0 && dtri2[1] < 0)) {
            std::swap(dtri2[0], dtri2[2]);
            std::swap(ptri2[0], ptri2[2]);
        } else if ((dtri2[0] > 0 && dtri2[2] > 0) || (dtri2[0] < 0 && dtri2[2] < 0)) {
            std::swap(dtri2[0], dtri2[1]);
            std::swap(ptri2[0], ptri2[1]);
        }
        
        // Line intervals for both triangles 
        t1_tri1 = ptri1[1] + (ptri1[0] - ptri1[1]) * dtri1[1] / (dtri1[1] - dtri1[0]);
        t2_tri1 = ptri1[2] + (ptri1[0] - ptri1[2]) * dtri1[2] / (dtri1[2] - dtri1[0]);
        t1_tri2 = ptri2[1] + (ptri2[0] - ptri2[1]) * dtri2[1] / (dtri2[1] - dtri2[0]);
        t2_tri2 = ptri2[2] + (ptri2[0] - ptri2[2]) * dtri2[2] / (dtri2[2] - dtri2[0]);
        
        // Check intervals intersection
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

    bool tri_tri3D(const triangle3D &tri1, const triangle3D &tri2) {
        return tri_tri3D(tri1.data(), tri2.data());
    }

    bool polyh_polyh(const polyhedron &polyhedron1, const polyhedron &polyhedron2) {
        for (const auto &tri1 : polyhedron1)
            for (const auto &tri2 : polyhedron2)
                if (tri_tri3D(tri1, tri2))
                    return true;
        return false;
    }

    // Check intersection between two polyherdons (using tri_tri3D function)
    //--------------------------------------------------------------------------------------------
    bool polyh_polyh(const Vector<3> (*_polyh1)[3], int _num_faces1, const Vector<3> (*_polyh2)[3], int _num_faces2)
    {
        for (int i = 0; i < _num_faces1; i++)
            for (int j = 0; j < _num_faces2; j++)
                if (tri_tri3D(_polyh1[i], _polyh2[j]))
                    return true;
        return false;
    }

    // Check whether the point lies inside convex polyhedron with positive-oriented faces
    //--------------------------------------------------------------------------------------------
    bool point_polyh(const Vector<3> & _point, const Vector<3>(*_polyh)[3], int _num_faces)
    {
        for (int i = 0; i < _num_faces; i++)
        {
            Vector<3> N = (_polyh[i][1] - _polyh[i][0]) ^ (_polyh[i][2] - _polyh[i][1]);
            double p = -(N * _polyh[i][0]);
            // Signed distances must be non-positive
            if (N * _point + p > 0)
                return false;
        }
        return true;
    }

    bool line_circle(const Vector<2> &center, double r, const Vector<2> &p1, const Vector<2> &p2)
    {
        double a = p2[0] - p1[0];
        double b = p2[1] - p1[1];
        double c = center[0] - p1[0];
        double d = center[1] - p1[1];
        if ((d*a - c*b) * (d*a - c*b) <= r*r * (a*a + b*b)) {
            if (c*c + d*d <= r*r) // first end
                return true;
            else if ((a - c) * (a - c) + (b - d) * (b - d) <= r*r) // second end
                return true;
            else if ((c*a + d*b >= 0) && (c*a + d*b <= a*a + b*b))
                return true;
            else
                return false;
        } else {
            return false;
        }
    }

    Vector<2> line_line(const Vector<2> &line1, double line1Angle, const Vector<2> &line2, double line2Angle)
    {
        double a = tan(line1Angle);
        double b = tan(line2Angle);
        double c = line1[1] - a * line1[0];
        double d = line2[1] - b * line2[0];

        double a_b = a - b;
        if (a_b == 0) // parallel lines
            return Vector<2>{{line1[0], line1[1]}};
        else
            return Vector<2>{{(d - c) / a_b, (a * d - b * c) / a_b}};
    }
}
