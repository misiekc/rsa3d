//
// Created by ciesla on 17.12.2025.
//

#ifndef RSA3D_ORIENTATIONSMODE_H
#define RSA3D_ORIENTATIONSMODE_H


#include "../ProgramMode.h"
#include "../ProgramArguments.h"

class OrientationsMode : public ProgramMode {
private:
    std::string packingFilename;

public:
    void initializeForArguments(const ProgramArguments &arguments) override;

    void printHelp(std::ostream &out, const ProgramArguments &arguments) override;

    void run() override;
};


#endif //RSA3D_ORIENTATIONSMODE_H
