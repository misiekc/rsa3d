//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_BCEXPANDMODE_H
#define RSA3D_BCEXPANDMODE_H


#include "../ProgramMode.h"
#include "../ProgramArguments.h"

class BCExpandMode : public ProgramMode {
private:
    std::string packingInFilename;
    std::string packingOutFilename;

public:
    void initializeForArguments(const ProgramArguments &arguments) override;
    void run() override;
};


#endif //RSA3D_BCEXPANDMODE_H
