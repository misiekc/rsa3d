//----------------------------------------------------------------------------
// Class representing convex polyhedron in R^3, derived after Shape class.
//----------------------------------------------------------------------------
// (C)PKua 2017
//----------------------------------------------------------------------------

#ifndef _CONVEH_POLYHEDRON_H
    #define _CONVEX_POLYHEDRON_H

#include <vector>
#include <map>
#include <memory>

#include "../Shape.h"
#include "../Matrix.h"
#include "../Vector.h"
#include "../OrientedFace.h"


class ConvexPolyhedron : public Shape
{
private:
    typedef std::shared_ptr<Vector<3>>          vptr;
    typedef std::shared_ptr<OrientedFace>       fptr;

    std::vector<vptr>       vertices;           // Shape vertices
    std::vector<fptr>       faces;              // All faces of polyhedron
    std::map<vptr, std::vector<fptr>>    vert_faces_map;     // Map of faces containing vertex

protected:

    double smallest_dist;       // The smallest possible distance from center to surface

public:
    ConvexPolyhedron(const std::vector<Vector<3>> & _vertices);
    ~ConvexPolyhedron();

    void makeFace(const std::vector<std::size_t> & _vertices);

	double getNeighbourListCellSize();
    double getVoxelSize();
    void translate(double* v);
    int overlap(BoundaryConditions *bc, Shape *s);
    double getVolume();
    int pointInside(BoundaryConditions *bc, double* da);
};

#endif // _CONVEX_POLYHEDRON_H
