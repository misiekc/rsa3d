//
// Created by PKua on 11.08.18.
//

#ifndef RSA3D_DISTRUBUTION_H
#define RSA3D_DISTRUBUTION_H

#include "../../rsa3d/utils/Utils.h"
#include "../../rsa3d/RND.h"

class Distribution {
public:
    virtual RSAVector randomVector(RND *rnd) const = 0;
    virtual std::string getDescription() const = 0;
};

#endif //RSA3D_DISTRUBUTION_H
