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
double          Cuboid::volume;
unsigned char   Cuboid::staticDimension;    
double          Cuboid::neighbourListCellSize;
double          Cuboid::voxelSize;

// Default constructor creating new Cuboid in (0, 0, 0) with size set in
// Cuboid::initClass
//----------------------------------------------------------------------------
Cuboid::Cuboid(unsigned char dimension, const Matrix & rotation) : Shape(dimension), rotation(rotation)
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
    for (unsigned char i = 0; i < staticDimension; i++) {
        args_stream >> size[i];
        if (!args_stream)           throw std::runtime_error("Cuboid::initClass: invalid or missing dimentions");
        else if (size[i] <= 0.0)    throw std::runtime_error("Cuboid::initClass: non-positive size: " + std::to_string(size[i]));
    }
        
    // Calculate static params
    volume = std::accumulate(size, size + staticDimension, 1.0, std::multiplies<double>());
    neighbourListCellSize = std::sqrt(
        std::accumulate(size, size + staticDimension, 0.0, 
            [](double sum, double x){
                return sum + x * x;
            }));    // diagonal lenght
    voxelSize = 0.0;
    
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

// Method creating (dynamically alocated) cuboid with random orientation
//----------------------------------------------------------------------------
Shape * Cuboid::create(RND *rnd)
{
    Cuboid * cuboid = new Cuboid(staticDimension, Matrix::rotation3D(
        rnd->nextValue() * 2 * M_PI,
        std::asin(rnd->nextValue() * 2 - 1),
        rnd->nextValue() * 2 * M_PI));
        
#ifdef CUBOID_DEBUG
    std::cout << "Creating Cuboid:" << std::endl;
    std::cout << cuboid->rotation;
#endif

    return cuboid;
}

double Cuboid::getNeighbourListCellSize()
{
    return Cuboid::neighbourListCellSize;
}

double Cuboid::getVoxelSize()
{
    return Cuboid::voxelSize;
}

int Cuboid::overlap(BoundaryConditions *bc, Shape *s)
{
    return false;
}

double Cuboid::getVolume()
{
    return Cuboid::volume;
}

int Cuboid::pointInside(BoundaryConditions *bc, double* da)
{
    return false;
}
