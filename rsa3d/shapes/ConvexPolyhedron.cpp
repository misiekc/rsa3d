//----------------------------------------------------------------------------
// Class representing convex polyhedron in R^3, derived after Shape class.
//----------------------------------------------------------------------------
// (C)PKua 2017
//----------------------------------------------------------------------------


#include <algorithm>
#include <iterator>
#include <iostream>


#include "ConvexPolyhedron.h"


// Constructs polyhedron from given vertices. After construction no faces
// are defined, so one have to use ConvexPolyhedron::makeFace
//----------------------------------------------------------------------------
ConvexPolyhedron::ConvexPolyhedron(const std::vector<Vector<3>> & _vertices) : Shape(3)
{
    this->vertices.reserve(_vertices.size());
    
    // Copy vectors to internal storage vector - shared_ptrs are user for GC
    std::transform(_vertices.begin(), _vertices.end(),
        std::back_inserter(this->vertices),
        [](const Vector<3> & v) { 
            std::shared_ptr<Vector<3>> p = std::make_shared<Vector<3>>(v);
            std::cout << p.get() << " " << &v << std::endl;
            return p;
        });
}


// Trivial destructor
//----------------------------------------------------------------------------
ConvexPolyhedron::~ConvexPolyhedron()
{
    // Let C++ and std::shared_ptr do GC
}


// Makes face from enumerated vetrices - numbers in vector correspond to
// vertices order in cunstructor, starting from 0
//----------------------------------------------------------------------------
void ConvexPolyhedron::makeFace(const std::vector<std::size_t> & _vertices)
{
    std::vector<vptr> vert_ptrs(_vertices.size());
    
    // Insert Vector<3> ptrs from shared_ptrs to vert_ptrs
    std::transform(_vertices.begin(), _vertices.end(),
        std::back_inserter(vert_ptrs),
        [&](std::size_t idx) {
            return vertices[idx];
        });
    
    // Add face to vector
    auto face = std::make_shared<OrientedFace>(vert_ptrs);
    this->faces.push_back(face);
    
    // Add face to vertices map
    for (auto v : vert_ptrs)
        this->vert_faces_map[v].push_back(face);
}


double ConvexPolyhedron::getNeighbourListCellSize()
{
    return 0;
}


double ConvexPolyhedron::getVoxelSize()
{
    return 0;
}


// Translates center and every vertex by vector v
//----------------------------------------------------------------------------
void ConvexPolyhedron::translate(double* v)
{
    Shape::translate(v);
    Vector<3> translation(v);
    for (auto v : this->vertices)
        *v += translation;
}


int ConvexPolyhedron::overlap(BoundaryConditions *bc, Shape *s)
{
    return 0;
}


// Returns volume of polyhedron
//----------------------------------------------------------------------------
double ConvexPolyhedron::getVolume()
{
    double vol = 0;
    Vector<3> center(this->position);
    // Minus, because signed distance is negative - center is inside polyhedron.
    // Divide polyhedron into pyramids with apices in center and bases on faces
    // and sum their volume
    for (auto f : this->faces)
        vol -= f->pointSignedDistance(center) * f->getSurface() / 3;
    return vol;
}


int ConvexPolyhedron::pointInside(BoundaryConditions *bc, double* da)
{
    Vector<3> point(da);
    for (auto f : this->faces)
        if (f->pointSignedDistance(point) > smallest_dist)
            return (int)false;
    return (int)true;
}

