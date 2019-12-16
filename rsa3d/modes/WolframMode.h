//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_WOLFRAMMODE_H
#define RSA3D_WOLFRAMMODE_H


#include "../ProgramMode.h"
#include "../ProgramArguments.h"

class WolframMode : public ProgramMode {
private:
    std::string packingFilename;
    bool isPeriodicImage{};

public:
    void initializeForArguments(const ProgramArguments &arguments) override;
    void printHelp(std::ostream &out, const ProgramArguments &arguments) override;
    void run() override;
};


#endif //RSA3D_WOLFRAMMODE_H
