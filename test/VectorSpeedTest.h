//
// Created by PKua on 04.03.18.
//

#ifndef RSA3D_VECTORSPEEDTEST_H
#define RSA3D_VECTORSPEEDTEST_H

#include "utility/Timer.h"
#include "../rsa3d/Vector.h"

namespace vec_speedtest {
    int main(int argc, char **argv);

    void printResult(const double *array1, Timer &timer);

    void printResult(const Vector<2> &vec1, Timer &timer);
}


#endif //RSA3D_VECTORSPEEDTEST_H
