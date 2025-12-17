//
// Created by user on 16.12.2025.
//

#include "OrientationAnalyzer.h"
#include "../shape/Shape.h"

void OrientationAnalyzer::analyzeOrientations(const Packing &packing) {
    for (const RSAShape *s: packing) {

        RSAOrientation orientation = s->getOrientation();
        Matrix<3,3> rotation = Matrix<3,3>::rotation(orientation[0], orientation[1], orientation[2]);
        Vector<3> direction = rotation*Vector<3>({1,0,0}).normalized();
        std::cout << direction[0] << ", " << direction[1] << ", " << direction[2] << std::endl;
    }
}
