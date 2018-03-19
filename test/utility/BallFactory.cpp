//--------------------------------------------------------------------------------------------
// Factory generating cuboids from ball of given radius
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "BallFactory.h"
#include "../../rsa3d/ShapeFactory.h"

typedef ShapePairFactory::ShapePair pair;
BallFactory * BallFactory::instance = nullptr;


// Private singleton constructor
//--------------------------------------------------------------------------------------------  
BallFactory::BallFactory()
{
    
}

// Returns singleton instance
//--------------------------------------------------------------------------------------------  
BallFactory * BallFactory::getInstance()
{
    if (instance == nullptr)
        instance = new BallFactory();
    return instance;
}


// Helper method. Creates random Cuboid based on objects parameters. Delegate cuboid creation
// to standard ShapeFactory
//--------------------------------------------------------------------------------------------    
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


// Sets sizes of half lengths of intervals describing box
//--------------------------------------------------------------------------------------------
void BallFactory::setRadius(double _radius)
{
    this->radius = _radius;
}


// Generates random pair of cuboids from the box of set size
//--------------------------------------------------------------------------------------------
pair BallFactory::generate()
{
    return pair(
            this->randomShape(),
        this->randomShape());
}


// Prints description of the factory on the standard output
//--------------------------------------------------------------------------------------------
std::string BallFactory::getDescription()
{
    std::stringstream stream;
    stream << "BallFactory of radius " << this->radius;
    return stream.str();
}


