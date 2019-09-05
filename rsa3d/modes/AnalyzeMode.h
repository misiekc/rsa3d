//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_ANALYZEMODE_H
#define RSA3D_ANALYZEMODE_H


#include "../ProgramMode.h"
#include "../ProgramArguments.h"

class AnalyzeMode : public ProgramMode {
private:
    std::string dirName;
    double corrRange{};
public:
    explicit AnalyzeMode(const ProgramArguments &arguments);

    void run() override;
};


#endif //RSA3D_ANALYZEMODE_H
