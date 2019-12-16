//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_DEBUGMODE_H
#define RSA3D_DEBUGMODE_H


#include "../ProgramMode.h"
#include "../ProgramArguments.h"

class DebugMode : public ProgramMode {
private:
    std::string packingGeneratorFilename;

public:
    void initializeForArguments(const ProgramArguments &arguments) override;
    void run() override;
};


#endif //RSA3D_DEBUGMODE_H
