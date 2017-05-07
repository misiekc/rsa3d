//----------------------------------------------------------------------------
// Class representing cuboid, derived after Shape class.
//----------------------------------------------------------------------------
// (C)PKua 2017
//----------------------------------------------------------------------------

#include "Cuboid.h"

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <numeric>
#include <functional>
#include <algorithm>
#include <iterator>

#define _USE_MATH_DEFINES
#include <cmath>


// Static atributes declaration
//----------------------------------------------------------------------------

double          *Cuboid::size;
double          *Cuboid::auxDoubleArray;
double          Cuboid::volume;
unsigned char   Cuboid::staticDimension;    
double          Cuboid::neighbourListCellSize;
double          Cuboid::voxelSize;

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


// Default constructor creating new Cuboid in (0, 0, 0) with size set in
// Cuboid::initClass
//----------------------------------------------------------------------------
Cuboid::Cuboid(const Matrix & orientation) : Shape(Cuboid::staticDimension), orientation(orientation)
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
    auxDoubleArray = new double[staticDimension];
    for (unsigned char i = 0; i < staticDimension; i++) {
        args_stream >> size[i];
        if (!args_stream)           throw std::runtime_error("Cuboid::initClass: invalid or missing dimentions");
        else if (size[i] <= 0.0)    throw std::runtime_error("Cuboid::initClass: non-positive size: " + std::to_string(size[i]));
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
    voxelSize = *std::min_element(size, size + staticDimension) / 2 / std::sqrt(staticDimension);
    
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
    Cuboid * cuboid = new Cuboid(Matrix::rotation3D(
        rnd->nextValue() * 2 * M_PI,
        std::asin(rnd->nextValue() * 2 - 1),
        rnd->nextValue() * 2 * M_PI));
        
#ifdef CUBOID_DEBUG
    std::cout << "Creating Cuboid:" << std::endl;
    std::cout << cuboid->orientation;
#endif

    return cuboid;
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

// Checks whether there is an overlap between this and *s
//----------------------------------------------------------------------------
int Cuboid::overlap(BoundaryConditions *bc, Shape *s)
{
    // Prepare matrices of translations for operations on shapes
    Cuboid *sCuboid = (Cuboid*)s;
    Matrix thisTranslation(this->dimension, 1, this->position);
    Matrix sTranslation(this->dimension, 1, sCuboid->position);
    sTranslation += Matrix(this->dimension, 1, bc->getTranslation(auxDoubleArray, this->position, sCuboid->position));
    Matrix backwards_rot = this->orientation.transpose();
    
    // Transfotm s coordinates to this coordinate system
    Matrix new_orientation = backwards_rot * sCuboid->orientation;
    sTranslation = backwards_rot * (sTranslation - thisTranslation);
    
    Matrix v_trans[V::SIZE];    // Calculated vertices coordinates
    
    // Check whether vertices of s lie in this. TO OPTIMIZE
    if (checkPoint( (v_trans[V::PPP] = new_orientation * Matrix(3, 1, { size[C::X] / 2,  size[C::Y] / 2,  size[C::Z] / 2}) + sTranslation) ) ||      
        checkPoint( (v_trans[V::NPP] = new_orientation * Matrix(3, 1, {-size[C::X] / 2,  size[C::Y] / 2,  size[C::Z] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::PNP] = new_orientation * Matrix(3, 1, { size[C::X] / 2, -size[C::Y] / 2,  size[C::Z] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::PPN] = new_orientation * Matrix(3, 1, { size[C::X] / 2,  size[C::Y] / 2, -size[C::Z] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::PNN] = new_orientation * Matrix(3, 1, { size[C::X] / 2, -size[C::Y] / 2, -size[C::Z] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::NPN] = new_orientation * Matrix(3, 1, {-size[C::X] / 2,  size[C::Y] / 2, -size[C::Z] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::NNP] = new_orientation * Matrix(3, 1, {-size[C::X] / 2, -size[C::Y] / 2,  size[C::Z] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::NNN] = new_orientation * Matrix(3, 1, {-size[C::X] / 2, -size[C::Y] / 2, -size[C::Z] / 2}) + sTranslation) ))
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
    return false;
}

// Checks whether given vertex (in this coordinates) lies in Cuboid
//----------------------------------------------------------------------------
bool Cuboid::checkPoint(const Matrix & vertex)
{
    for (unsigned char i = 0; i < this->dimension; i++)
        if (std::abs(vertex(i, 0)) > this->size[i] / 2)
            return false;
    return true;
}

// Checks whether a segment determined by point1 and point2 intersects with
// Cuboid
//----------------------------------------------------------------------------
bool Cuboid::checkSegment(const Matrix & point1, const Matrix & point2)
{
    double hsize_x = size[C::X] / 2;
    double hsize_y = size[C::Y] / 2;
    double hsize_z = size[C::Z] / 2;

    // Check intersections with all faces
    return checkSegmentFace( hsize_x, hsize_y, hsize_z, point1(C::X, 0), point1(C::Y, 0), point1(C::Z, 0), point2(C::X, 0), point2(C::Y, 0), point2(C::Z, 0)) ||
           checkSegmentFace(-hsize_x, hsize_y, hsize_z, point1(C::X, 0), point1(C::Y, 0), point1(C::Z, 0), point2(C::X, 0), point2(C::Y, 0), point2(C::Z, 0)) ||
           checkSegmentFace( hsize_y, hsize_z, hsize_x, point1(C::Y, 0), point1(C::Z, 0), point1(C::X, 0), point2(C::Y, 0), point2(C::Z, 0), point2(C::X, 0)) ||
           checkSegmentFace(-hsize_y, hsize_z, hsize_x, point1(C::Y, 0), point1(C::Z, 0), point1(C::X, 0), point2(C::Y, 0), point2(C::Z, 0), point2(C::X, 0)) ||
           checkSegmentFace( hsize_z, hsize_x, hsize_y, point1(C::Z, 0), point1(C::X, 0), point1(C::Y, 0), point2(C::Z, 0), point2(C::X, 0), point2(C::Y, 0)) ||
           checkSegmentFace(-hsize_z, hsize_x, hsize_y, point1(C::Z, 0), point1(C::X, 0), point1(C::Y, 0), point2(C::Z, 0), point2(C::X, 0), point2(C::Y, 0));          
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

// Checks whether the point with coordinates da lies inside Cuboid
//----------------------------------------------------------------------------
int Cuboid::pointInside(BoundaryConditions *bc, double* da)
{
    // Prepare matrices of translations for operations on the point
    Matrix cuboidTranslation(this->dimension, 1, this->position);
    Matrix pointTranslation(this->dimension, 1, da);
    
    // Transform point coordinates to Cuboid coordinate system
    cuboidTranslation += Matrix(this->dimension, 1, bc->getTranslation(auxDoubleArray, this->position, da));
    pointTranslation = this->orientation.transpose() * (pointTranslation - cuboidTranslation);
    
    // Check whether the point lies inside Cuboid
    return this->checkPoint(pointTranslation);
}
