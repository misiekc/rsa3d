//--------------------------------------------------------------------------------------------
// Class representing positive oriented face of convex 3D polyhedron
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include <stdexcept>
#include <iterator>
#include <algorithm>

#include "OrientedFace.h"


// Constructor taking vector of vertices (as R^3 Vectors). The caller is obligated to assure
// that veries are co-planar and lie in counter-clockwise order looking from "positive" side
//--------------------------------------------------------------------------------------------
OrientedFace::OrientedFace(const std::vector<Vector<3> *> & _vertices) : vertices(_vertices)
{
    if (vertices.size() < 3)
        throw std::runtime_error("face with less than 3 vertices");
}


// Constructor taking array of vertices (as R^3 Vectors). The caller is obligated to assure
// that veries are co-planar and lie in counter-clockwise order looking from "positive" side
//--------------------------------------------------------------------------------------------
OrientedFace::OrientedFace(Vector<3> ** _vertices, std::size_t _num)
{
    if (_num < 3)
        throw std::runtime_error("face with less than 3 vertices");
    
    this->vertices.reserve(_num);
    std::copy(_vertices, _vertices + _num, std::back_inserter(this->vertices));
}


// Returns vector orthogonal to face, positive oriented
//--------------------------------------------------------------------------------------------
Vector<3> OrientedFace::getOrthogonal() const
{
    return (*(this->vertices[1]) - *(this->vertices[0])) ^ (*(this->vertices[2]) - *(this->vertices[1]));
}


// Returns normal vector (normalized orthogonal)
//--------------------------------------------------------------------------------------------
Vector<3> OrientedFace::getNormal() const
{
    Vector<3> ort = this->getOrthogonal();
    return ort / ort.norm();
}


// Returns distance from point to face with corresponding sign
//--------------------------------------------------------------------------------------------
double OrientedFace::pointSignedDistance(const Vector<3> & _point) const
{
    Vector<3> N = this->getNormal();
    double d = -(N * *(this->vertices[0]));
    return N * _point + d;
}


// Returns negative number if point lies "under" the face, positive it "above" and zero
// if on
//--------------------------------------------------------------------------------------------
double OrientedFace::pointSign(const Vector<3> & _point) const
{
    Vector<3> N = this->getOrthogonal();
    double d = -(N * *(this->vertices[0]));
    return N * _point + d;
}


// Translates all vertices by given vector
//--------------------------------------------------------------------------------------------
void OrientedFace::translate(const Vector<3> & _translation)
{
    std::for_each(this->vertices.begin(), this->vertices.end(), [&](Vector<3> * v){ (*v) += _translation; });
}


// Returns vertices of this face (as vector of R^3 Vector)
//--------------------------------------------------------------------------------------------
std::vector<Vector<3> *> OrientedFace::getVertices() const
{
    return this->vertices;
}
