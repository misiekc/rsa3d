//
// Created by PKua on 11.08.18.
//

#include <iterator>
#include <algorithm>
#include "UniformBoxDistribution.h"

RSAVector UniformBoxDistribution::randomVector(RND *rnd) const {
    std::array<double, RSA_SPATIAL_DIMENSION> trans{};
    std::transform(this->halfsizes.begin(), this->halfsizes.end(), trans.begin(), [rnd](double hs) {
        return (rnd->nextValue() * 2 - 1) * hs;
    });
    return trans;
}

std::string UniformBoxDistribution::getDescription() const {
    std::ostringstream stream;
    stream << "UniformBoxDistribution of size ";
    std::copy(this->halfsizes.begin(), this->halfsizes.end(), std::ostream_iterator<double>(stream, " x "));
    return stream.str();
}

UniformBoxDistribution::UniformBoxDistribution(const std::array<double, RSA_SPATIAL_DIMENSION> &halfsizes)
        : halfsizes(halfsizes) {
    for (auto size : halfsizes)
        if (size < 0)
            throw std::runtime_error("one of halfsizes < 0");
}

UniformBoxDistribution::UniformBoxDistribution(double halfsize) {
    if (halfsize < 0)
        throw std::runtime_error("halfsize < 0");
    halfsizes.fill(halfsize);
}
