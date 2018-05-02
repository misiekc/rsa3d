//--------------------------------------------------------------------------------------------
// A package of functions for detecting intersections between primitives
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _INTERSECTION_H
    #define _INTERSECTION_H

#include "Matrix.h"
#include "Vector.h"
#include <vector>

namespace intersection
{
    using triangle3D = std::array<Vector<3>, 3>;
    using polyhedron = std::vector<triangle3D>;

    bool tri_tri3D(const triangle3D &tri1, const triangle3D &tri2);
    bool polyh_polyh(const polyhedron &polyhedron1, const polyhedron &polyhedron2);
    bool point_polyh(const Vector<3> & _point, const Vector<3> (*_polyh)[3], int _num_faces2);
    bool line_circle(const Vector<2> &center, double r, const Vector<2> &p1, const Vector<2> &p2);
    Vector<2> line_line(const Vector<2> &line1, double line1Angle, const Vector<2> &line2, double line2Angle);
}

#endif // _INTERSECTION_H
