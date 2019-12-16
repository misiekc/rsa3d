//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_PROGRAMMODE_H
#define RSA3D_PROGRAMMODE_H

#include "ProgramArguments.h"

class ProgramMode {
protected:
    Parameters params;

public:
    virtual ~ProgramMode() = default;

    virtual void initializeForArguments(const ProgramArguments &arguments) = 0;
    virtual void run() = 0;
    virtual void printHelp(std::ostream &out, const ProgramArguments &arguments) = 0;
};


#endif //RSA3D_PROGRAMMODE_H
