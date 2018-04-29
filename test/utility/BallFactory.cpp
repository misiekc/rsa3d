//--------------------------------------------------------------------------------------------
// Factory generating cuboids from ball of given radius
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include <sstream>

#include "BallFactory.h"
#include "../../rsa3d/ShapeFactory.h"

typedef ShapePairFactory::ShapePair pair;


Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> * BallFactory::randomShape()
{
    double trans[3];
    Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> * shape = ShapeFactory::createShape(&this->rnd);
    shape->no = this->no++;
    
    double cos_theta = 2 * rnd.nextValue() - 1;
    double sin_theta = sqrt(1 - cos_theta * cos_theta);
    double phi = rnd.nextValue() * 2 * M_PI;
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);
    double radius = pow(rnd.nextValue(), 1. / 3) * this->radius;
    
    trans[0] = radius * sin_theta * cos_phi;
    trans[1] = radius * sin_theta * sin_phi;
    trans[2] = radius * cos_theta;
    shape->translate(trans);
    return shape;
}

void BallFactory::setRadius(double _radius) { this->radius = _radius; }

pair BallFactory::generate() { return {this->randomShape(), this->randomShape()}; }

std::string BallFactory::getDescription() const
{
    std::stringstream stream;
    stream << "BallFactory of radius " << this->radius;
    return stream.str();
}


