//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_PACKINGTESTMODE_H
#define RSA3D_PACKINGTESTMODE_H


#include "../ProgramMode.h"
#include "../ProgramArguments.h"

class PackingTestMode : public ProgramMode {
private:
    std::string packingFilename;
    double maxTime{};

public:
    void initializeForArguments(const ProgramArguments &arguments) override;
    void printHelp(std::ostream &out, const ProgramArguments &arguments) override;
    void run() override;
};


#endif //RSA3D_PACKINGTESTMODE_H
