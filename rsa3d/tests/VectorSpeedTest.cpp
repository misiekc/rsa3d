//
// Created by PKua on 04.03.18.
//

#include <iostream>
#include "VectorSpeedTest.h"
#include "../RND.h"

#define REPEATS 100000000

namespace
{
    /* Functions for testing passing and returning. All in two copies. */
    void addArray1(double *result, double *arr1, double *arr2) {
        for (std::size_t j = 0; j < 2; j++)
            result[j] = arr1[j] + arr2[j];
    }

    void addArray2(double *result, double *arr1, double *arr2) {
        for (std::size_t j = 0; j < 2; j++)
            result[j] = arr1[j] + arr2[j];
    }

    Vector<2> addVector1(const Vector<2> &vec1, const Vector<2> &vec2) {
        return vec1 + vec2;
    }

    Vector<2> addVector2(const Vector<2> &vec1, const Vector<2> &vec2) {
        return vec1 + vec2;
    }

    void addVectorReference1(Vector<2> &result, const Vector<2> &vec1, const Vector<2> &vec2) {
        result = vec1 + vec2;
    }

    void addVectorReference2(Vector<2> &result, const Vector<2> &vec1, const Vector<2> &vec2) {
        result = vec1 + vec2;
    }

    void reset(double (&array1)[2], double (&array2)[2], double (&array3)[2],
               Vector<2> &vec1, Vector<2> &vec2, Vector<2> &vec3) {
        array1[0] = 0.1; array1[1] = 0.1;
        array2[0] = 0.04; array2[1] = 0.04;
        array3[0] = 0; array3[1] = 0;

        vec1 = Vector<2>(array1);
        vec2 = Vector<2>(array2);
        vec3 = Vector<2>();
    }
}

using arr_func_type = decltype(&addArray1);
using vec_func_type = decltype(&addVector1);
using vec_ref_func_type = decltype(&addVectorReference1);


int vec_speedtest::main(int argc, char **argv) {
    Timer timer;
    RND rnd;

    // Prevent inlining. Compiler would even inline function pointed by a function pointer if it can be
    // determined in compile-time. Wow
    arr_func_type arr_func = rnd.nextValue() > 0.5 ? addArray1 : addArray2;
    vec_func_type vec_func = rnd.nextValue() > 0.5 ? addVector1 : addVector2;
    vec_ref_func_type vec_ref_func = rnd.nextValue() > 0.5 ? addVectorReference1 : addVectorReference2;

    double array1[2];
    double array2[2];
    double array3[2];
    Vector<2> vec1;
    Vector<2> vec2;
    Vector<2> vec3;

    std::cout << ">> Warm up..." << std::endl;
    reset(array1, array2, array3, vec1, vec2, vec3);
    timer.start();
    for (std::size_t i = 0; i < 4000000000; i++)
        for (std::size_t j = 0; j < 2; j++)
            array1[j] += array2[j];
    timer.stop();
    printResult(array1, timer);
    std::cout << std::endl;

    std::cout << ">> Addition:" << std::endl;
    reset(array1, array2, array3, vec1, vec2, vec3);
    timer.start();
    for (std::size_t i = 0; i < REPEATS; i++)
        for (std::size_t j = 0; j < 2; j++)
            array1[j] += array2[j];
    timer.stop();
    printResult(array1, timer);

    timer.start();
    for (std::size_t i = 0; i < REPEATS; i++)
        vec1 += vec2;
    timer.stop();
    printResult(vec1, timer);
    std::cout << std::endl;

    std::cout << ">> Addition by function:" << std::endl;
    reset(array1, array2, array3, vec1, vec2, vec3);
    double arrayTmp[2];
    timer.start();
    for (std::size_t i = 0; i < REPEATS; i++) {
        arr_func(arrayTmp, array1, array2);
        for (std::size_t j = 0; j < 2; j++)
            array3[j] += arrayTmp[j];
    }
    timer.stop();
    printResult(array3, timer);

    Vector<2> vecTmp;
    timer.start();
    for (std::size_t i = 0; i < REPEATS; i++) {
        vec_ref_func(vecTmp, vec1, vec2);
        vec3 += vecTmp;
    }
    timer.stop();
    std::cout << "(reference) ";
    printResult(vec3, timer);

    reset(array1, array2, array3, vec1, vec2, vec3);
    timer.start();
    for (std::size_t i = 0; i < REPEATS; i++)
        vec3 += vec_func(vec1, vec2);
    timer.stop();
    std::cout << "(returned) ";
    printResult(vec3, timer);
    std::cout << std::endl;

    return EXIT_SUCCESS;
}

inline void vec_speedtest::printResult(const Vector<2> &vec1, Timer &timer) {
    std::cout << "vector " << vec1 << ": " << timer.count<std::chrono::milliseconds>() << " ms" << std::endl;
}

inline void vec_speedtest::printResult(const double *array1, Timer &timer) {
    std::cout << "array {" << array1[0] << ", " << array1[1] << "}: " << timer.count<std::chrono::milliseconds>() << " ms" << std::endl;
}
