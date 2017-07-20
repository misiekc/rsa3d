//----------------------------------------------------------------------------
// Class representing convex polyhedron in R^3, derived after Shape class.
//----------------------------------------------------------------------------
// (C)PKua 2017
//----------------------------------------------------------------------------


#include <algorithm>
#include <iterator>
#include <limits>
#include <iostream>
#include <sstream>


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
            return new Vector<3>(v);
        });
}


// Cunstructs polyhedron from given cuboid
//----------------------------------------------------------------------------
ConvexPolyhedron::ConvexPolyhedron(Cuboid * _cube) : Shape(3)
{
    Vector<3> pos(_cube->getPosition()); 
    Matrix<3, 3> orientation = _cube->getOrientation();
    double size[3];
    Cuboid::getSize(size);
    
    // Calculate verties
    this->vertices.reserve(8);
    this->vertices.push_back(new Vector<3>(pos + orientation * Vector<3>{{ size[0] / 2,  size[1] / 2,  size[2] / 2}}));
    this->vertices.push_back(new Vector<3>(pos + orientation * Vector<3>{{-size[0] / 2,  size[1] / 2,  size[2] / 2}}));
    this->vertices.push_back(new Vector<3>(pos + orientation * Vector<3>{{ size[0] / 2, -size[1] / 2,  size[2] / 2}}));
    this->vertices.push_back(new Vector<3>(pos + orientation * Vector<3>{{ size[0] / 2,  size[1] / 2, -size[2] / 2}}));
    this->vertices.push_back(new Vector<3>(pos + orientation * Vector<3>{{ size[0] / 2, -size[1] / 2, -size[2] / 2}}));
    this->vertices.push_back(new Vector<3>(pos + orientation * Vector<3>{{-size[0] / 2,  size[1] / 2, -size[2] / 2}}));
    this->vertices.push_back(new Vector<3>(pos + orientation * Vector<3>{{-size[0] / 2, -size[1] / 2,  size[2] / 2}}));
    this->vertices.push_back(new Vector<3>(pos + orientation * Vector<3>{{-size[0] / 2, -size[1] / 2, -size[2] / 2}}));
    
    // Make faces
    this->faces.reserve(6);
    this->makeFace(std::vector<size_t>{0, 3, 5, 1});
    this->makeFace(std::vector<size_t>{0, 2, 4, 3});
    this->makeFace(std::vector<size_t>{2, 6, 7, 4});
    this->makeFace(std::vector<size_t>{6, 1, 5, 7});
    this->makeFace(std::vector<size_t>{0, 1, 6, 2});
    this->makeFace(std::vector<size_t>{4, 7, 5, 3});
}


// Destructor
//----------------------------------------------------------------------------
ConvexPolyhedron::~ConvexPolyhedron()
{
    for (auto f : this->faces)
        delete f;
    for (auto v : this->vertices)
        delete v;
}


// Makes face from enumerated vetrices - numbers in vector correspond to
// vertices order in cunstructor, starting from 0
//----------------------------------------------------------------------------
void ConvexPolyhedron::makeFace(const std::vector<std::size_t> & _vertices)
{
    std::vector<vptr> vert_ptrs;
    vert_ptrs.reserve(_vertices.size());
    
    // Insert Vector<3> ptrs from shared_ptrs to vert_ptrs
    std::transform(_vertices.begin(), _vertices.end(),
        std::back_inserter(vert_ptrs),
        [&](std::size_t idx) {
            return vertices[idx];
        });
    
    // Add face to vector
    fptr face = new OrientedFace(vert_ptrs);
    this->faces.push_back(face);
    
    // Add face to vertices map
    /*for (auto v : vert_ptrs)
        this->vert_faces_map[v].push_back(face);*/
}


double ConvexPolyhedron::getNeighbourListCellSize()
{
    return 0;
}


double ConvexPolyhedron::getVoxelSize()
{
    return 0;
}


// Translates center and every vertex by vector v (double array)
//----------------------------------------------------------------------------
void ConvexPolyhedron::translate(double* _v)
{
    Shape::translate(_v);
    Vector<3> translation(_v);
    for (auto vert : this->vertices)
        *vert += translation;
}


// Translates center and every vertex by vector v (Vector<3>)
//----------------------------------------------------------------------------
void ConvexPolyhedron::translate(const Vector<3> & _v)
{
    double v_arr[3];
    for (size_t i = 0; i < 3; i++)
        v_arr[i] = _v[i];
    Shape::translate(v_arr);
    for (auto vert : this->vertices)
        *vert += _v;
}

    
// Checks whether polyhedron overlaps with another one
//----------------------------------------------------------------------------
int ConvexPolyhedron::overlap(BoundaryConditions *bc, Shape *s)
{
    ConvexPolyhedron * sPolyh = (ConvexPolyhedron *)s;
    double trans_arr[3];
    Vector<3> translation(bc->getTranslation(trans_arr, this->position, sPolyh->position));
    
    // Translate shape s for PBC
    sPolyh->translate(translation);
    
    // Check all possible separating axes - normals ...
    for (auto f : this->faces)
        if (!checkAxis (f->getOrthogonal(), sPolyh))
            goto ret_false;
    for (auto f : sPolyh->faces)
        if (!checkAxis (f->getOrthogonal(), sPolyh))
            goto ret_false;
    // ... and normals cross products
    for (auto f1 : this->faces)
        for (auto f2 : sPolyh->faces)
            if (!checkAxis (f1->getOrthogonal() ^ f2->getOrthogonal(), sPolyh))
                goto ret_false;
   
    // Translate back before return
    sPolyh->translate(-translation);
    return true;
    
    ret_false:
    sPolyh->translate(-translation);
    return false;
}


// Checks whether this and _second projections on axis _axis overlap. If so,
// returns true
//----------------------------------------------------------------------------
bool ConvexPolyhedron::checkAxis(const Vector<3> & _axis, const ConvexPolyhedron * _second) const
{
    interval this_int = this->getProjection(_axis, this);
    interval second_int = this->getProjection(_axis, _second);
    
    return std::min(this_int.second, second_int.second) >= std::max(this_int.first, second_int.first);
}


// Projects polyhedron _polyh on axis _axis and returns interval given by
// the projection
//----------------------------------------------------------------------------
ConvexPolyhedron::interval ConvexPolyhedron::getProjection(const Vector<3> & _axis, const ConvexPolyhedron * _polyh)
{
    interval proj_int = {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    
    // Find enpoints of polyhedron projection (multiplied by unknown but const for _axit factor)
    double proj;
    for (auto v : _polyh->vertices)
    {
        proj = *v * _axis;
        if (proj < proj_int.first)
            proj_int.first = proj;
        if (proj > proj_int.second)
            proj_int.second = proj; 
    }
    return proj_int;
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


// Check weather point lies in excluded volume
//----------------------------------------------------------------------------
int ConvexPolyhedron::pointInside(BoundaryConditions *bc, double* da)
{
    Vector<3> point(da);
    for (auto f : this->faces)
        if (f->pointSignedDistance(point) > smallest_dist)
            return (int)false;
    return (int)true;
}


// Returns Wolfram Language representation
//----------------------------------------------------------------------------
std::string ConvexPolyhedron::toWolfram() const
{
    std::stringstream stream;
    
    stream << "convex_poly" << this->no << " = Polygon[{" << std::endl;
    for (std::size_t i = 0; i < this->faces.size() - 1; i++) {
        this->wolframFaceToStream(this->faces[i], stream);
        stream << "," << std::endl;
    }
    this->wolframFaceToStream(this->faces.back(), stream);
    stream << "}];" << std::endl;
    
    return stream.str();
}


// Prints face to stringstream. Helper for ConvexPolyhedron::toWolfram
//----------------------------------------------------------------------------
void ConvexPolyhedron::wolframFaceToStream(fptr _face, std::stringstream & _stream) const
{
    _stream << "    {";
    std::vector<vptr> vs = _face->getVertices();
    for (std::size_t i = 0; i < vs.size() - 1; i++)
        _stream << *(vs[i]) << ", ";
    _stream << *(vs.back());
    _stream << "}";
}


