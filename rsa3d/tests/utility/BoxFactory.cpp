//--------------------------------------------------------------------------------------------
// Factory generating cuboids from box of given size
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "BoxFactory.h"
#include "../../ShapeFactory.h"
#include <sstream>


typedef ShapePairFactory::ShapePair pair;
BoxFactory * BoxFactory::instance = nullptr;


// Private singleton constructor
//--------------------------------------------------------------------------------------------  
BoxFactory::BoxFactory()
{
    
}

// Returns singleton instance
//--------------------------------------------------------------------------------------------  
BoxFactory * BoxFactory::getInstance()
{
    if (instance == nullptr)
        instance = new BoxFactory();
    return instance;
}


// Helper method. Creates random Cuboid based on objects parameters. Delegate cuboid creation
// to standard ShapeFactory
//--------------------------------------------------------------------------------------------    
Shape<RSA_DIMENSION> * BoxFactory::randomShape()
{
    double trans[3];
    Shape<RSA_DIMENSION> * shape = ShapeFactory::createShape(&this->rnd);
    shape->no = this->no++;
    trans[0] = (rnd.nextValue() * 2 - 1) * this->halfsizeX;
    trans[1] = (rnd.nextValue() * 2 - 1) * this->halfsizeY;
    trans[2] = (rnd.nextValue() * 2 - 1) * this->halfsizeZ;
    shape->translate(trans);
    return shape;
}


// Sets sizes of half lengths of intervals describing box
//--------------------------------------------------------------------------------------------
void BoxFactory::setBoxSize(double _halfsize_x, double _halfsize_y, double _halfsize_z)
{
    this->halfsizeX = _halfsize_x;
    this->halfsizeY = _halfsize_y;
    this->halfsizeZ = _halfsize_z;
}


// Generates random pair of cuboids from the box of set size
//--------------------------------------------------------------------------------------------
pair BoxFactory::generate()
{
    return pair(
            this->randomShape(),
        this->randomShape());
}


// Prints description of the factory on the standard output
//--------------------------------------------------------------------------------------------
std::string BoxFactory::getDescription()
{
    std::stringstream stream;
    stream << "BoxFactory of size " << this->halfsizeX * 2 << " x " << this->halfsizeY * 2 
        << " x " << this->halfsizeZ * 2;
    return stream.str();
}


