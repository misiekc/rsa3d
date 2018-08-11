//
// Created by PKua on 11.08.18.
//

#ifndef RSA3D_BOXUNIFORMDISTRIBUTION_H
#define RSA3D_BOXUNIFORMDISTRIBUTION_H


#include "Distrubution.h"

class UniformBoxDistribution : public Distribution {
private:
    std::array<double, RSA_SPATIAL_DIMENSION> halfsizes{};

public:
    UniformBoxDistribution() : UniformBoxDistribution(0) {}
    explicit UniformBoxDistribution(const std::array<double, RSA_SPATIAL_DIMENSION> &halfsizes);
    explicit UniformBoxDistribution(double halfsize);

    RSAVector randomVector(RND *rnd) const override;
    std::string getDescription() const override;
};


#endif //RSA3D_BOXUNIFORMDISTRIBUTION_H
