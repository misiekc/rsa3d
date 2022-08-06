//
// Created by pkua on 06.11.2021.
//

#ifndef RSA3D_MYTESTMODE_H
#define RSA3D_MYTESTMODE_H

#include "../ProgramMode.h"

class MyTestMode : public ProgramMode {
public:
    void run() override;
    void printHelp(std::ostream &out, const ProgramArguments &arguments) override;
};


#endif //RSA3D_MYTESTMODE_H
