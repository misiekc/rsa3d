/*
 * Utils.c
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */


#ifndef UTILS_C_
#define UTILS_C_

#include "Vector.h"

#ifdef _OPENMP
    #define __OMP_STRINGIFY__(x) #x

    #define _OMP_PARALLEL_FOR   _Pragma("omp parallel for")
    #define _OMP_ATOMIC         _Pragma("omp atomic")
    #define _OMP_CRITICAL(x)    _Pragma(__OMP_STRINGIFY__(omp critical(x)))
    #define _OMP_MAXTHREADS     omp_get_max_threads()
    #define _OMP_THREAD_ID      omp_get_thread_num()
#else
    #define _OMP_PARALLEL_FOR
    #define _OMP_ATOMIC
    #define _OMP_CRITICAL(x)
    #define _OMP_MAXTHREADS     1
    #define _OMP_THREAD_ID      0
#endif

using RSAVector = Vector<RSA_SPATIAL_DIMENSION>;

template <std::size_t N>
using Orientation = std::array<double, N>;

using RSAOrientation = Orientation<RSA_ANGULAR_DIMENSION>;

bool increment(int* in, int inlength, int max);
int position2i(const double* da, int dalength, double size, double dx, int n);
void i2position(double* da, int dalength, int index, double dx, int n);
void coordinates(int* result, const double* da, int dalength, double size, double dx, int n);
int neighbour2i(int* coordinates, int* neighbour, int clength, int offset, int n);

// trim from start
std::string &ltrim(std::string &s);
// trim from end
std::string &rtrim(std::string &s);
// trim from both ends
std::string &trim(std::string &s);
// replaces all occurences of search in source by replace
std::string replaceAll(std::string source, const std::string& search, const std::string& replace);

int lastIndexOf(const char * s, char target);

void die(const std::string & reason);

double getAngleToOrigin(const Vector<2> &point);
void rotate2D(double* point, double alpha);

#endif /* UTILS_C_ */
