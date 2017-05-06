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
//----------------------------------------------------------------------------
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
    if (checkPoint( (v_trans[V::PPP] = new_orientation * Matrix(3, 1, { size[0] / 2,  size[1] / 2,  size[2] / 2}) + sTranslation) ) ||      
        checkPoint( (v_trans[V::NPP] = new_orientation * Matrix(3, 1, {-size[0] / 2,  size[1] / 2,  size[2] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::PNP] = new_orientation * Matrix(3, 1, { size[0] / 2, -size[1] / 2,  size[2] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::PPN] = new_orientation * Matrix(3, 1, { size[0] / 2,  size[1] / 2, -size[2] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::PNN] = new_orientation * Matrix(3, 1, { size[0] / 2, -size[1] / 2, -size[2] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::NPN] = new_orientation * Matrix(3, 1, {-size[0] / 2,  size[1] / 2, -size[2] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::NNP] = new_orientation * Matrix(3, 1, {-size[0] / 2, -size[1] / 2,  size[2] / 2}) + sTranslation) ) ||
        checkPoint( (v_trans[V::NNN] = new_orientation * Matrix(3, 1, {-size[0] / 2, -size[1] / 2, -size[2] / 2}) + sTranslation) ))
    {    
        return true;    
    }
    
    // Check whether edges of s lie in this
    
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
