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

using Triangle3D = std::array<Vector<3>, 3>;
using Face3D = std::vector<Vector<3>>;
using PolyhedronTriangulation = std::vector<Triangle3D>;
using PolyhedronFaces = std::vector<Face3D>;

namespace collision
{
    bool tri3D_tri3D(const Triangle3D &tri1, const Triangle3D &tri2);
    bool polyh_polyh(const PolyhedronTriangulation &polyhedron1, const PolyhedronTriangulation &polyhedron2);
    bool point_polyh(const Vector<3> & _point, const Vector<3> (*_polyh)[3], int _num_faces2);
    bool line_circle(const Vector<2> &center, double r, const Vector<2> &p1, const Vector<2> &p2);
}

namespace intersection
{
    Vector<2> line_line(const Vector<2> &line1, double line1Angle, const Vector<2> &line2, double line2Angle);
}

#endif // _INTERSECTION_H
