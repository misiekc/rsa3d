//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_PROGRAMMODE_H
#define RSA3D_PROGRAMMODE_H

#include "Parameters.h"

class ProgramMode {
protected:
    Parameters params;

public:
    explicit ProgramMode(Parameters params) : params{std::move(params)} { }

    virtual void run() = 0;
};


#endif //RSA3D_PROGRAMMODE_H
