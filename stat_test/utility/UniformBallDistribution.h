//
// Created by PKua on 11.08.18.
//

#ifndef RSA3D_BALLUNIFORMDISTRIBUTION_H
#define RSA3D_BALLUNIFORMDISTRIBUTION_H


#include "Distrubution.h"
#include "../../rsa3d/RND.h"

class UniformBallDistribution : public Distribution {
private:
    double radius;

    double randomGaussian(RND *rnd) const;

public:
    UniformBallDistribution() : UniformBallDistribution(0) {}
    explicit UniformBallDistribution(double radius);

    RSAVector randomVector(RND *rnd) const override;
    std::string getDescription() const override;

    void setRadius(double radius);
};


#endif //RSA3D_BALLUNIFORMDISTRIBUTION_H
