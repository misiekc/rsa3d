//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_DATFILEGENERATIONMODE_H
#define RSA3D_DATFILEGENERATIONMODE_H


#include "../ProgramMode.h"
#include "../ProgramArguments.h"

class DatFileGenerationMode : public ProgramMode {
private:
    std::string dirName;

    double getMedian(std::vector<double> &data);

public:
    explicit DatFileGenerationMode(const ProgramArguments &arguments);

    void run() override;
};


#endif //RSA3D_DATFILEGENERATIONMODE_H
