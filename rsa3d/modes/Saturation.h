//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_SATURATION_H
#define RSA3D_SATURATION_H


#include "../Packing.h"
#include "../Parameters.h"
#include "../ProgramMode.h"

class Saturation : public ProgramMode {

private:
    std::string packingFilename;

	unsigned long getSeedOrigin() const;

public:
    void initializeForArguments(const ProgramArguments &arguments) final;
    void run() final;
    void printHelp(std::ostream &out, const ProgramArguments &arguments) override;
};


#endif //RSA3D_SATURTATION_H
