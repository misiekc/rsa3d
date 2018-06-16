//--------------------------------------------------------------------------------------------
// Factory generating cuboids from box of given size
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include <sstream>
#include <iterator>
#include "BoxFactory.h"
#include "../../rsa3d/ShapeFactory.h"
#include "ShapeGenerators.h"


typedef RSAShapePairFactory::ShapePair pair;

BoxFactory::BoxFactory(double halfsize) {
    this->halfsize.fill(halfsize);
}

std::unique_ptr<RSAShape> BoxFactory::randomShape() {
    RSAShape *shape = ShapeFactory::createShape(&this->rnd);
    shape->no = this->no++;

    std::array<double, RSA_SPATIAL_DIMENSION> trans{};
    std::transform(this->halfsize.begin(), this->halfsize.end(), trans.begin(), [this](double hs) {
        return (rnd.nextValue() * 2 - 1) * hs;
    });

    return generate_randomly_oriented_shape(Vector<RSA_SPATIAL_DIMENSION>(trans), &rnd);
}

void BoxFactory::setBoxSize(const std::array<double, RSA_SPATIAL_DIMENSION> &halfsize) {
    this->halfsize = halfsize;
}

pair BoxFactory::generate() { return {this->randomShape(), this->randomShape()}; }

std::string BoxFactory::getDescription() const {
    std::ostringstream stream;
    stream << "BoxFactory of size ";
    std::copy(this->halfsize.begin(), this->halfsize.end(), std::ostream_iterator<double>(stream, " x "));
    return stream.str();
}


