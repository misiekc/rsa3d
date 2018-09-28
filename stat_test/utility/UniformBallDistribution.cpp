//
// Created by PKua on 11.08.18.
//

#include <algorithm>
#include "UniformBallDistribution.h"

UniformBallDistribution::UniformBallDistribution(double radius) : radius(radius) {
    if (radius < 0.0)  throw std::runtime_error("radius < 0");
}

RSAVector UniformBallDistribution::randomVector(RND *rnd) const {
    std::array<double, RSA_SPATIAL_DIMENSION> posArray{};
    std::for_each(posArray.begin(), posArray.end(), [rnd, this](double &elem){
        elem = this->randomGaussian(rnd);
    });

    double radius = pow(rnd->nextValue(), 1. / RSA_SPATIAL_DIMENSION) * this->radius;
    Vector<RSA_SPATIAL_DIMENSION> pos(posArray);
    return pos.normalized() * radius;
}

double UniformBallDistribution::randomGaussian(RND *rnd) const {
    return std::sqrt(-2 * std::log(rnd->nextValue())) * std::cos(2 * M_PI * rnd->nextValue());
}

std::string UniformBallDistribution::getDescription() const {
    return "UniformBallDistribution of radius " + std::to_string(this->radius);
}

void UniformBallDistribution::setRadius(double radius) {
    if (radius < 0)    throw std::runtime_error("radius < 0");
    this->radius = radius;
}
