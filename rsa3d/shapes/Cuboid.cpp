//----------------------------------------------------------------------------
// Class representing cuboid, derived after Shape class.
//----------------------------------------------------------------------------
// (C)PKua 2017
//----------------------------------------------------------------------------

#include "Cuboid.h"
#include "../Vector.h"
#include "../Intersection.h"
#include "ConvexPolyhedron.h"

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <numeric>
#include <functional>
#include <algorithm>
#include <iterator>
#include <limits>

#define _USE_MATH_DEFINES
#include <cmath>


// Define this to disable edge-intersection fix in Cuboid::OverlapStrategy::MINE
//#define DISABLE_OVERLAP_FIX

namespace 
{
    // Helper arrays
    Vector<3>       cuboid1_tris[12][3];
    Vector<3>       cuboid2_tris[12][3];
}


// Static atributes declaration
//----------------------------------------------------------------------------

double          *Cuboid::size;
double          *Cuboid::auxDoubleArray;
double          Cuboid::volume;
unsigned char   Cuboid::staticDimension;    
double          Cuboid::neighbourListCellSize;
double          Cuboid::voxelSize;
double          Cuboid::minDimension;

// Vertex recognition helper. P states positive, N - negative. First position
// corresponds to positive/negative X, second for Y, etc.
enum V {
    PPP = 0,
    NPP,
    PNP,
    PPN,
    PNN,
    NPN,
    NNP,
    NNN,
    SIZE
};

// Coord recognition helper
enum C {
    X = 0,
    Y,
    Z
};

inline bool checkSegmentFace(double plane_pos, double bound1, double bound2, 
                             double p1p, double p1b1, double p1b2,
                             double p2p, double p2b1, double p2b2);


Cuboid::OverlapStrategy Cuboid::strategy = Cuboid::OverlapStrategy::MINE;


// Default constructor creating new Cuboid in (0, 0, 0) with size set in
// Cuboid::initClass
//----------------------------------------------------------------------------
Cuboid::Cuboid(const Matrix<3, 3> & orientation) : Shape(Cuboid::staticDimension), orientation(orientation)
{
    
}

// Destructor (does nothing)
//----------------------------------------------------------------------------
Cuboid::~Cuboid()
{

}

// Static method for class initialization
//----------------------------------------------------------------------------
// args - string arguments. Format: <size x> <size y> <size z>
//----------------------------------------------------------------------------
void Cuboid::initClass(const std::string &args)
{
    std::stringstream args_stream (args);
    int temp;
    
    // Parse parameters
    args_stream >> temp;
    staticDimension = (unsigned char)temp;
    if (!args_stream)               throw std::runtime_error("Cuboid::initClass: invalid arguments. Usage: <dimension> <size 1> ... <size dimension>");
    else if (staticDimension == 0)  throw std::runtime_error("Cuboid::initClass: 0 dimension");
        
    size = new double[staticDimension];
    auxDoubleArray = new double[staticDimension];
    for (unsigned char i = 0; i < staticDimension; i++) {
        args_stream >> size[i];
        if (!args_stream)           throw std::runtime_error("Cuboid::initClass: invalid or missing dimensions");
        else if (size[i] <= 0.0)    throw std::runtime_error("Cuboid::initClass: non-positive size: " + std::to_string(size[i]));
    }

    // renormailze sizes to obtain unit volume
    double v = std::accumulate(size, size + staticDimension, 1.0, std::multiplies<double>());
    double factor = 1.0/pow(v, 1.0/staticDimension);
    for (unsigned char i = 0; i < staticDimension; i++) {
        size[i] *= factor;
    }
        
    // Calculate static params
    volume = std::accumulate(size, size + staticDimension, 1.0, std::multiplies<double>());
    
    // Neighbour list cell size - diagonal length. It satisfies conditions
    // at the least favourable case - Cuboids with centers near opposite edges
    // of a cell and with diagonals lying on the line connecting their centers
    neighbourListCellSize = std::sqrt(
        std::accumulate(size, size + staticDimension, 0.0, 
            [](double sum, double x){
                return sum + x * x;
            }));
            
    // Voxel size. It satisfies conditions at the least favourable case
    // - center of Cuboid lies in the corner of a voxel and voxel's diagonal
    // is parallel to the smallest Cuboid edge
    minDimension = *std::min_element(size, size + staticDimension) / 2;
    voxelSize = minDimension / std::sqrt(staticDimension);
    
#ifdef CUBOID_DEBUG
    std::cout << "Initializing Cuboid class" << std::endl;
    std::cout << "dimensions           : ";
    // Implode dimensions with comma delimiter an write to standard output
    std::copy(size, size + staticDimension, std::ostream_iterator<double>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "volume               : " << volume << std::endl;
    std::cout << "cell size (diagonal) : " << neighbourListCellSize << std::endl;
#endif
}

// Method creating (dynamically alocated) cuboid with random orientation.
// Used by ShapeFactory for shape generating
//----------------------------------------------------------------------------
Shape * Cuboid::create(RND *rnd)
{
    Cuboid * cuboid = new Cuboid(Matrix<3, 3>::rotation(
        rnd->nextValue() * 2 * M_PI,
        std::asin(rnd->nextValue() * 2 - 1),
        rnd->nextValue() * 2 * M_PI));
        
#ifdef CUBOID_DEBUG
    std::cout << "Creating Cuboid:" << std::endl;
    std::cout << cuboid->orientation;
#endif

    return cuboid;
}


// Sets which overlap algorithm to use
//----------------------------------------------------------------------------
void Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy _strategy)
{
    strategy = _strategy;
}


// Returns overlap algorithm used
//----------------------------------------------------------------------------
Cuboid::OverlapStrategy Cuboid::getOverlapStrategy()
{
    return Cuboid::strategy;
}


// Helper method. Obtains and saves triangles from cuboid's faces
//--------------------------------------------------------------------------------------------
void Cuboid::obtainTris(Vector<3> (&arr)[12][3], const Vector<3> & translation)
{
    Vector<3> pos(this->getPosition()); 
    Matrix<3, 3> orientation = this->getOrientation();   
    Vector<3> vert[] = {
        pos + translation + orientation * Vector<3>{{ size[C::X] / 2,  size[C::Y] / 2,  size[C::Z] / 2}},
        pos + translation + orientation * Vector<3>{{-size[C::X] / 2,  size[C::Y] / 2,  size[C::Z] / 2}},
        pos + translation + orientation * Vector<3>{{ size[C::X] / 2, -size[C::Y] / 2,  size[C::Z] / 2}},
        pos + translation + orientation * Vector<3>{{ size[C::X] / 2,  size[C::Y] / 2, -size[C::Z] / 2}},
        pos + translation + orientation * Vector<3>{{ size[C::X] / 2, -size[C::Y] / 2, -size[C::Z] / 2}},
        pos + translation + orientation * Vector<3>{{-size[C::X] / 2,  size[C::Y] / 2, -size[C::Z] / 2}},
        pos + translation + orientation * Vector<3>{{-size[C::X] / 2, -size[C::Y] / 2,  size[C::Z] / 2}},
        pos + translation + orientation * Vector<3>{{-size[C::X] / 2, -size[C::Y] / 2, -size[C::Z] / 2}} 
    };
    
    arr[0][0] = vert[0];
    arr[0][1] = vert[1];
    arr[0][2] = vert[2];

    arr[1][0] = vert[2];
    arr[1][1] = vert[1];
    arr[1][2] = vert[6];

    arr[2][0] = vert[0];
    arr[2][1] = vert[2];
    arr[2][2] = vert[3];

    arr[3][0] = vert[3];
    arr[3][1] = vert[2];
    arr[3][2] = vert[4];

    arr[4][0] = vert[7];
    arr[4][1] = vert[2];
    arr[4][2] = vert[6];

    arr[5][0] = vert[7];
    arr[5][1] = vert[4];
    arr[5][2] = vert[2];

    arr[6][0] = vert[1];
    arr[6][1] = vert[0];
    arr[6][2] = vert[3];

    arr[7][0] = vert[5];
    arr[7][1] = vert[1];
    arr[7][2] = vert[3];

    arr[8][0] = vert[7];
    arr[8][1] = vert[1];
    arr[8][2] = vert[5];

    arr[9][0] = vert[7];
    arr[9][1] = vert[6];
    arr[9][2] = vert[1];

    arr[10][0] = vert[7];
    arr[10][1] = vert[5];
    arr[10][2] = vert[3];

    arr[11][0] = vert[7];
    arr[11][1] = vert[3];
    arr[11][2] = vert[4];
}


// Returns neighbour list cell size determined during class initialization
//----------------------------------------------------------------------------
double Cuboid::getNeighbourListCellSize()
{
    return Cuboid::neighbourListCellSize;
}

// Returns initial voxel size determined during class initialization
//----------------------------------------------------------------------------
double Cuboid::getVoxelSize()
{
    return Cuboid::voxelSize;
}

// Checks whether there is an overlap between this and *s using chosen
// stategy
//----------------------------------------------------------------------------
int Cuboid::overlap(BoundaryConditions *bc, Shape *s)
{
    switch (strategy) {
        case TRI_TRI:
            return overlapTri(bc, s);
        case SAT:
            return overlapSAT(bc, s);
        default:
            return overlapMine(bc, s);
    }
}

// My overlap algorithm
//----------------------------------------------------------------------------
int Cuboid::overlapMine(BoundaryConditions *bc, Shape *s)
{    
    // Prepare matrices of translations for operations on shapes
    Cuboid *sCuboid = (Cuboid*)s;
    Vector<3> thisTranslation(this->position);
    Vector<3> sTranslation(sCuboid->position);
    sTranslation += Vector<3>(bc->getTranslation(auxDoubleArray, this->position, sCuboid->position));
    Matrix<3, 3> backwards_rot = this->orientation.transpose();
    
    // Transform s coordinates to this coordinate system
    Matrix<3, 3> new_orientation = backwards_rot * sCuboid->orientation;
    sTranslation = backwards_rot * (sTranslation - thisTranslation);
    
    Vector<3> v_trans[V::SIZE];    // Calculated vertices coordinates
    Vector<3> v_trans_bis[V::SIZE];    // Calculated vertices coordinates in swapped cuboid order
    
    // Check whether vertices of s lie in this. TO OPTIMIZE
    if (checkPoint( (v_trans[V::PPP] = new_orientation * Vector<3>{{ size[C::X] / 2,  size[C::Y] / 2,  size[C::Z] / 2}} + sTranslation) ) ||      
        checkPoint( (v_trans[V::NPP] = new_orientation * Vector<3>{{-size[C::X] / 2,  size[C::Y] / 2,  size[C::Z] / 2}} + sTranslation) ) ||
        checkPoint( (v_trans[V::PNP] = new_orientation * Vector<3>{{ size[C::X] / 2, -size[C::Y] / 2,  size[C::Z] / 2}} + sTranslation) ) ||
        checkPoint( (v_trans[V::PPN] = new_orientation * Vector<3>{{ size[C::X] / 2,  size[C::Y] / 2, -size[C::Z] / 2}} + sTranslation) ) ||
        checkPoint( (v_trans[V::PNN] = new_orientation * Vector<3>{{ size[C::X] / 2, -size[C::Y] / 2, -size[C::Z] / 2}} + sTranslation) ) ||
        checkPoint( (v_trans[V::NPN] = new_orientation * Vector<3>{{-size[C::X] / 2,  size[C::Y] / 2, -size[C::Z] / 2}} + sTranslation) ) ||
        checkPoint( (v_trans[V::NNP] = new_orientation * Vector<3>{{-size[C::X] / 2, -size[C::Y] / 2,  size[C::Z] / 2}} + sTranslation) ) ||
        checkPoint( (v_trans[V::NNN] = new_orientation * Vector<3>{{-size[C::X] / 2, -size[C::Y] / 2, -size[C::Z] / 2}} + sTranslation) ))
    {    
        return true;    
    }
    
    // Transform this coordinates to s coordinate system
    new_orientation = new_orientation.transpose();
    sTranslation = -(new_orientation * sTranslation);
    
    // Check whether vertices of this lie in s
    if (sCuboid->checkPoint( (v_trans_bis[V::PPP] = new_orientation * Vector<3>{{ size[C::X] / 2,  size[C::Y] / 2,  size[C::Z] / 2}} + sTranslation) ) ||      
        sCuboid->checkPoint( (v_trans_bis[V::NPP] = new_orientation * Vector<3>{{-size[C::X] / 2,  size[C::Y] / 2,  size[C::Z] / 2}} + sTranslation) ) ||
        sCuboid->checkPoint( (v_trans_bis[V::PNP] = new_orientation * Vector<3>{{ size[C::X] / 2, -size[C::Y] / 2,  size[C::Z] / 2}} + sTranslation) ) ||
        sCuboid->checkPoint( (v_trans_bis[V::PPN] = new_orientation * Vector<3>{{ size[C::X] / 2,  size[C::Y] / 2, -size[C::Z] / 2}} + sTranslation) ) ||
        sCuboid->checkPoint( (v_trans_bis[V::PNN] = new_orientation * Vector<3>{{ size[C::X] / 2, -size[C::Y] / 2, -size[C::Z] / 2}} + sTranslation) ) ||
        sCuboid->checkPoint( (v_trans_bis[V::NPN] = new_orientation * Vector<3>{{-size[C::X] / 2,  size[C::Y] / 2, -size[C::Z] / 2}} + sTranslation) ) ||
        sCuboid->checkPoint( (v_trans_bis[V::NNP] = new_orientation * Vector<3>{{-size[C::X] / 2, -size[C::Y] / 2,  size[C::Z] / 2}} + sTranslation) ) ||
        sCuboid->checkPoint( (v_trans_bis[V::NNN] = new_orientation * Vector<3>{{-size[C::X] / 2, -size[C::Y] / 2, -size[C::Z] / 2}} + sTranslation) ))
    {
        return true;    
    }
    
    // Check whether edges of s lie in this. TO OPTIMIZE
    if (checkSegment(v_trans[V::PPP], v_trans[V::PPN]) || checkSegment(v_trans[V::PPN], v_trans[V::PNN]) ||
        checkSegment(v_trans[V::PNN], v_trans[V::PNP]) || checkSegment(v_trans[V::PNP], v_trans[V::PPP]) ||
        checkSegment(v_trans[V::NNN], v_trans[V::NNP]) || checkSegment(v_trans[V::NNP], v_trans[V::NPP]) ||
        checkSegment(v_trans[V::NPP], v_trans[V::NPN]) || checkSegment(v_trans[V::NPN], v_trans[V::NNN]) ||
        checkSegment(v_trans[V::PPP], v_trans[V::NPP]) || checkSegment(v_trans[V::PPN], v_trans[V::NPN]) ||
        checkSegment(v_trans[V::PNN], v_trans[V::NNN]) || checkSegment(v_trans[V::PNP], v_trans[V::NNP]))
    {
        return true;
    }
    
#ifndef DISABLE_OVERLAP_FIX
    // Check whether edges of this lie in s. TO OPTIMIZE
    if (sCuboid->checkSegment(v_trans_bis[V::PPP], v_trans_bis[V::PPN]) || sCuboid->checkSegment(v_trans_bis[V::PPN], v_trans_bis[V::PNN]) ||
        sCuboid->checkSegment(v_trans_bis[V::PNN], v_trans_bis[V::PNP]) || sCuboid->checkSegment(v_trans_bis[V::PNP], v_trans_bis[V::PPP]) ||
        sCuboid->checkSegment(v_trans_bis[V::NNN], v_trans_bis[V::NNP]) || sCuboid->checkSegment(v_trans_bis[V::NNP], v_trans_bis[V::NPP]) ||
        sCuboid->checkSegment(v_trans_bis[V::NPP], v_trans_bis[V::NPN]) || sCuboid->checkSegment(v_trans_bis[V::NPN], v_trans_bis[V::NNN]) ||
        sCuboid->checkSegment(v_trans_bis[V::PPP], v_trans_bis[V::NPP]) || sCuboid->checkSegment(v_trans_bis[V::PPN], v_trans_bis[V::NPN]) ||
        sCuboid->checkSegment(v_trans_bis[V::PNN], v_trans_bis[V::NNN]) || sCuboid->checkSegment(v_trans_bis[V::PNP], v_trans_bis[V::NNP]))
    {
        return true;
    }
#endif
    
    return false;
}

// Triangle overlap algorithm
//----------------------------------------------------------------------------
int Cuboid::overlapTri(BoundaryConditions *bc, Shape *s)
{    
    Cuboid * second = (Cuboid *)s;
    double trans_arr[3];
    Vector<3> translation(bc->getTranslation(trans_arr, this->position, second->position));

    this->obtainTris(cuboid1_tris, Vector<3>());
    second->obtainTris(cuboid2_tris, translation);
    
    return intersection::polyh_polyh(cuboid1_tris, 12, cuboid2_tris, 12);
}


// SAT overlap algorithm
//----------------------------------------------------------------------------
int Cuboid::overlapSAT(BoundaryConditions *bc, Shape *s)
{
    Cuboid * second = (Cuboid *)s;
    double trans_arr[3];
    Vector<3> translation(bc->getTranslation(trans_arr, this->position, second->position));
    Vector<3> pos1(this->position); 
    Vector<3> pos2(second->position);
    Matrix<3, 3> orientation1 = this->orientation;
    Matrix<3, 3> orientation2 = second->orientation;
    Vector<3> vertices1[8];
    Vector<3> vertices2[8];
    double size[3];
    
    Cuboid::getSize(size);
    
    // Calculate axes orthogonal to separating plane
    Vector<3> axes1[] = {
        orientation1 * Vector<3>{{1, 0, 0}},
        orientation1 * Vector<3>{{0, 1, 0}},
        orientation1 * Vector<3>{{0, 0, 1}}
    };
    Vector<3> axes2[] = {
        orientation2 * Vector<3>{{1, 0, 0}},
        orientation2 * Vector<3>{{0, 1, 0}},
        orientation2 * Vector<3>{{0, 0, 1}}
    };
    
    
    // Calculate verties
    vertices1[0] = Vector<3>(pos1 + orientation1 * Vector<3>{{ size[0] / 2,  size[1] / 2,  size[2] / 2}});
    vertices1[1] = Vector<3>(pos1 + orientation1 * Vector<3>{{-size[0] / 2,  size[1] / 2,  size[2] / 2}});
    vertices1[2] = Vector<3>(pos1 + orientation1 * Vector<3>{{ size[0] / 2, -size[1] / 2,  size[2] / 2}});
    vertices1[3] = Vector<3>(pos1 + orientation1 * Vector<3>{{ size[0] / 2,  size[1] / 2, -size[2] / 2}});
    vertices1[4] = Vector<3>(pos1 + orientation1 * Vector<3>{{ size[0] / 2, -size[1] / 2, -size[2] / 2}});
    vertices1[5] = Vector<3>(pos1 + orientation1 * Vector<3>{{-size[0] / 2,  size[1] / 2, -size[2] / 2}});
    vertices1[6] = Vector<3>(pos1 + orientation1 * Vector<3>{{-size[0] / 2, -size[1] / 2,  size[2] / 2}});
    vertices1[7] = Vector<3>(pos1 + orientation1 * Vector<3>{{-size[0] / 2, -size[1] / 2, -size[2] / 2}});
    
    vertices2[0] = Vector<3>(pos2 + translation + orientation2 * Vector<3>{{ size[0] / 2,  size[1] / 2,  size[2] / 2}});
    vertices2[1] = Vector<3>(pos2 + translation + orientation2 * Vector<3>{{-size[0] / 2,  size[1] / 2,  size[2] / 2}});
    vertices2[2] = Vector<3>(pos2 + translation + orientation2 * Vector<3>{{ size[0] / 2, -size[1] / 2,  size[2] / 2}});
    vertices2[3] = Vector<3>(pos2 + translation + orientation2 * Vector<3>{{ size[0] / 2,  size[1] / 2, -size[2] / 2}});
    vertices2[4] = Vector<3>(pos2 + translation + orientation2 * Vector<3>{{ size[0] / 2, -size[1] / 2, -size[2] / 2}});
    vertices2[5] = Vector<3>(pos2 + translation + orientation2 * Vector<3>{{-size[0] / 2,  size[1] / 2, -size[2] / 2}});
    vertices2[6] = Vector<3>(pos2 + translation + orientation2 * Vector<3>{{-size[0] / 2, -size[1] / 2,  size[2] / 2}});
    vertices2[7] = Vector<3>(pos2 + translation + orientation2 * Vector<3>{{-size[0] / 2, -size[1] / 2, -size[2] / 2}});
    
    // Check all possible separating axes - edge lines ...
    for (int i = 0; i < 3; i++)
        if (!this->checkSeparatingAxis (axes1[i], vertices1, vertices2))
            return false;
    for (int i = 0; i < 3; i++)
        if (!this->checkSeparatingAxis (axes2[i], vertices1, vertices2))
            return false;
    // ... and their cross products
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (!this->checkSeparatingAxis (axes1[i] ^ axes2[j], vertices1, vertices2))
                return false;
    return true;
}

// Checks whether this and _second projections on axis _axis overlap. If so,
// returns true
//----------------------------------------------------------------------------
bool Cuboid::checkSeparatingAxis(const Vector<3> & _axis, Vector<3> * _vert1, Vector<3> * _vert2) const
{
    interval this_int = this->getProjection(_axis, _vert1);
    interval second_int = this->getProjection(_axis, _vert2);
    
    return std::min(this_int.second, second_int.second) >= std::max(this_int.first, second_int.first);
}


// Projects polyhedron _polyh on axis _axis and returns interval given by
// the projection
//----------------------------------------------------------------------------
Cuboid::interval Cuboid::getProjection(const Vector<3> & _axis, Vector<3> * _vert) const
{
    interval proj_int = {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    
    // Find enpoints of polyhedron projection (multiplied by unknown but const for _axit factor)
    double proj;
    for (std::size_t i = 0; i < 8; i++)
    {
        proj = _vert[i] * _axis;
        if (proj < proj_int.first)
            proj_int.first = proj;
        if (proj > proj_int.second)
            proj_int.second = proj; 
    }
    return proj_int;
}


// Checks whether given vertex (in this coordinates) lies in Cuboid
//----------------------------------------------------------------------------
bool Cuboid::checkPoint(const Vector<3> & vertex)
{
    for (unsigned char i = 0; i < this->dimension; i++)
        if (std::abs(vertex[i]) > this->size[i] / 2)
            return false;
    return true;
}

// Checks whether a segment determined by point1 and point2 intersects with
// Cuboid
//----------------------------------------------------------------------------
bool Cuboid::checkSegment(const Vector<3> & point1, const Vector<3> & point2)
{
    double hsize_x = size[C::X] / 2;
    double hsize_y = size[C::Y] / 2;
    double hsize_z = size[C::Z] / 2;

    // Check intersections with all faces
    return checkSegmentFace( hsize_x, hsize_y, hsize_z, point1[C::X], point1[C::Y], point1[C::Z], point2[C::X], point2[C::Y], point2[C::Z]) ||
           checkSegmentFace(-hsize_x, hsize_y, hsize_z, point1[C::X], point1[C::Y], point1[C::Z], point2[C::X], point2[C::Y], point2[C::Z]) ||
           checkSegmentFace( hsize_y, hsize_z, hsize_x, point1[C::Y], point1[C::Z], point1[C::X], point2[C::Y], point2[C::Z], point2[C::X]) ||
           checkSegmentFace(-hsize_y, hsize_z, hsize_x, point1[C::Y], point1[C::Z], point1[C::X], point2[C::Y], point2[C::Z], point2[C::X]) ||
           checkSegmentFace( hsize_z, hsize_x, hsize_y, point1[C::Z], point1[C::X], point1[C::Y], point2[C::Z], point2[C::X], point2[C::Y]) ||
           checkSegmentFace(-hsize_z, hsize_x, hsize_y, point1[C::Z], point1[C::X], point1[C::Y], point2[C::Z], point2[C::X], point2[C::Y]);          
}

// Checks whether given segment intersects with axis-oriented rectangular face
// with middle on the perpendicular axis (OP)
//----------------------------------------------------------------------------
// plane_pos - face-determined plane coordinate on OP
// bound1 - positive face edge coordinate on the 1st axis perpendicular to OP (OB1)
// bound2 - positive face edge coordinate on the 2nd axis perpendicular to OP (OB2)
// p1p, p1b1, p1b2 - 1st segment point coordinates respectively on OP, OB1, OB2
// p2p, p2b1, p2b2 - 2nd segment point coordinates respectively on OP, OB1, OB2
//----------------------------------------------------------------------------
inline bool checkSegmentFace(double plane_pos, double bound1, double bound2, 
                             double p1p, double p1b1, double p1b2,
                             double p2p, double p2b1, double p2b2)
{
    // Check weather a plane determined by face lies between segment points
    if ((plane_pos > p1p && plane_pos > p2p) || (plane_pos < p1p && plane_pos < p2p))
        return false;
    
    // Intersect plane's face with segment's line
    double b1i = ((p1b1 - p2b1) * plane_pos + p1p * p2b1 - p2p * p1b1) / (p1p - p2p);
    double b2i = ((p1b2 - p2b2) * plane_pos + p1p * p2b2 - p2p * p1b2) / (p1p - p2p);
    
    // Check whether intersection point lies on face
    if (std::abs(b1i) > bound1 || std::abs(b2i) > bound2)
        return false;
    else
        return true;
}


// Returns volume of the Cuboid determined during class initialization
//----------------------------------------------------------------------------
double Cuboid::getVolume()
{
    return Cuboid::volume;
}

// Checks whether the point with coordinates da lies inside excluded volume.
// It is an interior of a set of point which distanse from
//----------------------------------------------------------------------------
int Cuboid::pointInside(BoundaryConditions *bc, double* da)
{
    // Prepare matrices of translations for operations on the point
    Vector<3> cuboidTranslation(this->position);
    Vector<3> pointTranslation(da);
    
    // Transform point coordinates to Cuboid coordinate system
    cuboidTranslation += Vector<3>(bc->getTranslation(auxDoubleArray, this->position, da));
    pointTranslation = this->orientation.transpose() * (pointTranslation - cuboidTranslation);
    
    // Save absolute values of point coords
    for (unsigned char i = 0; i < 3; i++)
        auxDoubleArray[i] = std::abs(pointTranslation[i]);
    
    // Map which coords lie in this and which lie in this + minDimension
    bool liesInSmaller[3];
    bool liesInBigger[3];
    for (unsigned char i = 0; i < 3; i++) {
        liesInSmaller[i] = (auxDoubleArray[i] <= this->size[i] / 2);
        liesInBigger[i] = (auxDoubleArray[i] <= this->size[i] / 2 + minDimension);
    }
    
    // Check optimistic cases - "lies in this" and "doesn't lie in this + minDimension" 
    if (!liesInBigger[0] || !liesInBigger[1] || !liesInBigger[2])       return false;
    if (liesInSmaller[0] && liesInSmaller[1] && liesInSmaller[2])       return true;
    
    // Check "pushed rectangles"
    if ((liesInSmaller[0] && liesInSmaller[1]) || (liesInSmaller[1] && liesInSmaller[2]) || (liesInSmaller[2] && liesInSmaller[0]))    
        return true;
        
    // Check cylinders on edges
    if (liesInSmaller[0] && std::pow(auxDoubleArray[1] - this->size[1] / 2, 2) + std::pow(auxDoubleArray[2] - this->size[2] / 2, 2) <= std::pow(minDimension, 2))
        return true;
    if (liesInSmaller[1] && std::pow(auxDoubleArray[2] - this->size[2] / 2, 2) + std::pow(auxDoubleArray[0] - this->size[0] / 2, 2) <= std::pow(minDimension, 2))
        return true;
    if (liesInSmaller[2] && std::pow(auxDoubleArray[0] - this->size[0] / 2, 2) + std::pow(auxDoubleArray[1] - this->size[1] / 2, 2) <= std::pow(minDimension, 2))
        return true;
    
    // Check spheres in vertices
    if (std::pow(auxDoubleArray[0] - this->size[0] / 2, 2) + std::pow(auxDoubleArray[1] - this->size[1] / 2, 2) +
        std::pow(auxDoubleArray[2] - this->size[2] / 2, 2) <= std::pow(minDimension, 2))
        return true;
        
    return false;
}

// Returns Cuboid orientation
//----------------------------------------------------------------------------
Matrix<3, 3> Cuboid::getOrientation() const
{
    return this->orientation;
}

// Returns Cuboid size
//----------------------------------------------------------------------------
double * Cuboid::getSize(double * arr)
{
    std::copy(size, size + staticDimension, arr);
    return arr;
}

std::string Cuboid::toPovray() const
{
	std::string s = "";
	if (this->dimension==2){
		double factor = 0.5;
		s += "  polygon {4, ";

		s += "< " + std::to_string(this->position[0] - factor*size[0]) + ", ";
		s += 		std::to_string(this->position[1] - factor*size[1]) + ", 0.1>, ";

		s += "< " + std::to_string(this->position[0] - factor*size[0]) + ", ";
		s += 		std::to_string(this->position[1] + factor*size[1]) + ", 0.1>, ";

		s += "< " + std::to_string(this->position[0] + factor*size[0]) + ", ";
		s += 		std::to_string(this->position[1] + factor*size[1]) + ", 0.1>, ";

		s += "< " + std::to_string(this->position[0] + factor*size[0]) + ", ";
		s += 		std::to_string(this->position[1] - factor*size[1]) + ", 0.1> ";

		s += "\n    texture { pigment { color Red } }\n  }\n";
	}
	else if (this->dimension==3){
		s += "  box { < ";
		for(unsigned char i=0; i<this->dimension; i++){
			s += std::to_string(-this->size[i]/2);
			if (i<this->dimension-1)
				s+= ", ";
		}
		s += ">, <";
		for(unsigned char i=0; i<this->dimension; i++){
			s += std::to_string(this->size[i]/2);
			if (i<this->dimension-1)
				s+= ", ";
		}
		s += ">\n";
		s += "    matrix < \n    ";
		for (unsigned char i=0; i<this->dimension; i++){
			for (unsigned char j=0; j<this->dimension; j++){
				s += std::to_string(this->orientation(i, j));
				if (j<this->dimension-1)
					s+= ", ";
			}
			s+= ",\n    ";
		}

		for(unsigned char i=0; i<this->dimension; i++){
			s += std::to_string(this->position[i]);
			if (i<this->dimension-1)
				s+= ", ";
		}
		s += "\n    >\n";
		s += "    texture { pigment { color Red } }\n  }\n";
	}
	return s;
}

// Returns Wolfram Language for drawing shape
//----------------------------------------------------------------------------
std::string Cuboid::toWolfram() const
{
    if (staticDimension != 3)
        throw std::runtime_error("Cuboid::toWolfram supports only 3D Cuboids");
    
    std::stringstream out;
    out << "cube" << this->no << " = " << std::endl;
    out << "    GeometricTransformation[" << std::endl;
    out << "        Cuboid[{" << (-this->size[0] / 2) << ", " << (-this->size[1] / 2) << ", " << (-this->size[2] / 2) << "}, ";
    out << "{" << (this->size[0] / 2) << ", " << (this->size[1] / 2) << ", " << (this->size[2] / 2) << "}]," << std::endl;
    out << "        AffineTransform[" << std::endl;
    out << "            {{{" << this->orientation(0, 0) << ", " << this->orientation(0, 1) << ", " << this->orientation(0, 2) << "}," << std::endl;
    out << "            {" << this->orientation(1, 0) << ", " << this->orientation(1, 1) << ", " << this->orientation(1, 2) << "}," << std::endl;
    out << "            {" << this->orientation(2, 0) << ", " << this->orientation(2, 1) << ", " << this->orientation(2, 2) << "}}," << std::endl;
    out << "            {" << this->position[0] << ", " << this->position[1] << ", " << this->position[2] << "}}]];";
    return out.str();
}

void Cuboid::store(std::ostream &f) const
{
	Shape::store(f);
	double d;
	for (unsigned char i=0; i<this->dimension; i++){
		for (unsigned char j=0; j<this->dimension; j++){
			d = this->orientation(i, j);
			f.write((char *)(&d), sizeof(double));
		}
	}
}

void Cuboid::restore(std::istream &f)
{
	Shape::restore(f);
	double d;
	for (unsigned char i=0; i<this->dimension; i++){
		for (unsigned char j=0; j<this->dimension; j++){
			f.read((char *)&d, sizeof(double));
			this->orientation(i, j) = d;
		}
	}
}

