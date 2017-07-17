//--------------------------------------------------------------------------------------------
// A package of functions for detecting intersections between primitives
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _INTERSECTION_H
    #define _INTERSECTION_H

#include "Matrix.h"
#include "Vector.h"

namespace intersection
{
    bool tri_tri3D(const Vector<3> * _tri1, const Vector<3> * _tri2);
    bool polyh_polyh(const Vector<3> (*_polyh1)[3], int _num_faces1, const Vector<3> (*_polyh2)[3], int _num_faces2);
    bool point_polyh(const Vector<3> & _point, const Vector<3> (*_polyh)[3], int _num_faces2);
}

#endif // _INTERSECTION_H
