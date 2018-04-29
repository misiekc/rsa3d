//--------------------------------------------------------------------------------------------
// Factory generating cuboids from box of given size
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include <sstream>
#include "BoxFactory.h"
#include "../../rsa3d/ShapeFactory.h"


typedef ShapePairFactory::ShapePair pair;

Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> * BoxFactory::randomShape()
{
    double trans[3];
    Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> * shape = ShapeFactory::createShape(&this->rnd);
    shape->no = this->no++;
    trans[0] = (rnd.nextValue() * 2 - 1) * this->halfsizeX;
    trans[1] = (rnd.nextValue() * 2 - 1) * this->halfsizeY;
    trans[2] = (rnd.nextValue() * 2 - 1) * this->halfsizeZ;
    shape->translate(trans);
    return shape;
}

void BoxFactory::setBoxSize(double _halfsize_x, double _halfsize_y, double _halfsize_z)
{
    this->halfsizeX = _halfsize_x;
    this->halfsizeY = _halfsize_y;
    this->halfsizeZ = _halfsize_z;
}

pair BoxFactory::generate() { return {this->randomShape(), this->randomShape()}; }

std::string BoxFactory::getDescription() const
{
    std::stringstream stream;
    stream << "BoxFactory of size " << this->halfsizeX * 2 << " x " << this->halfsizeY * 2 
        << " x " << this->halfsizeZ * 2;
    return stream.str();
}


