//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_EXCLUSIONZONESMODE_H
#define RSA3D_EXCLUSIONZONESMODE_H


#include "../ProgramMode.h"
#include "../ProgramArguments.h"

class ExclusionZonesMode : public ProgramMode {
private:
    std::string packingFilename;
    std::string outputFilename;

public:
    void initializeForArguments(const ProgramArguments &arguments) override;
    void run() override;
};


#endif //RSA3D_EXCLUSIONZONESMODE_H
