//
// Created by Michal Ciesla on 16.11.2022.
//

#ifndef RSA3D_DOMAINSMODE_H
#define RSA3D_DOMAINSMODE_H


#include "../ProgramMode.h"
#include "../ProgramArguments.h"

class DomainsMode : public ProgramMode {
private:
    std::string dirFile;

public:
    void initializeForArguments (const ProgramArguments &arguments) override;
    void printHelp(std::ostream &out, const ProgramArguments &arguments) override;
    void run() override;
};


#endif //RSA3D_DOMAINSMODE_H
