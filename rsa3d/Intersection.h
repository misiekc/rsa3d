//--------------------------------------------------------------------------------------------
// A package of functions for detecting intersections between primitives
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "Matrix.h"
#include "Vector.h"

namespace intersection
{
    bool tri_tri3D(const Vector * _tri1, const Vector * _tri2);
    bool polyh_polyh(const Vector (*_polyh1)[3], int _num_faces1, const Vector (*_polyh2)[3], int _num_faces2);
    bool point_polyh(const Vector & _point, const Vector(*_polyh)[3], int _num_faces2);
}
