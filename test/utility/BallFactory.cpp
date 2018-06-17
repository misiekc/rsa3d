//--------------------------------------------------------------------------------------------
// Factory generating cuboids from ball of given radius
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include <sstream>

#include "BallFactory.h"
#include "../../rsa3d/shape/ShapeFactory.h"
#include "ShapeGenerators.h"

typedef RSAShapePairFactory::ShapePair pair;

std::unique_ptr<RSAShape> BallFactory::randomShape() {
    RSAShape *shape = ShapeFactory::createShape(&this->rnd);
    shape->no = this->no++;

    std::array<double, RSA_SPATIAL_DIMENSION> posArray{};
    std::for_each(posArray.begin(), posArray.end(), [this](double &elem){
       elem = this->randomGaussian();
    });

    double radius = pow(rnd.nextValue(), 1. / RSA_SPATIAL_DIMENSION) * this->radius;
    Vector<RSA_SPATIAL_DIMENSION> pos(posArray);
    pos = pos / pos.norm() * radius;

    return generate_randomly_oriented_shape(pos, &rnd);
}

void BallFactory::setRadius(double _radius) { this->radius = _radius; }

pair BallFactory::generate() { return {this->randomShape(), this->randomShape()}; }

std::string BallFactory::getDescription() const {
    std::stringstream stream;
    stream << "BallFactory of radius " << this->radius;
    return stream.str();
}

double BallFactory::randomGaussian() {
    return std::sqrt(-2 * std::log(rnd.nextValue())) * std::cos(2 * M_PI * rnd.nextValue());
}


