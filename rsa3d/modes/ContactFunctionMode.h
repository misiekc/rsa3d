//
// Created by pkua on 05.09.2019.
//

#ifndef RSA3D_CONTACTFUNCTIONMODE_H
#define RSA3D_CONTACTFUNCTIONMODE_H


#include "../ProgramMode.h"
#include "../ProgramArguments.h"
#include "../shape/Shape.h"
#include "../boundary_conditions/FreeBC.h"

class ContactFunctionMode : public ProgramMode {
private:
	RSAFreeBC bc;
	int stepAngle{};
    void absoluteTranslate(RSAShape *s, RSAVector *v);
    double calculate(RSAShape *s1, RSAShape *s2);
    void calculate2D(std::string filename);

public:
    void initializeForArguments(const ProgramArguments &arguments) override;
    void printHelp(std::ostream &out, const ProgramArguments &arguments) override;
    void run() override;
};


#endif //RSA3D_CONTACTFUNCTIONMODE_H
