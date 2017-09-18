/*
 * Utils.c
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */


#ifndef UTILS_C_
#define UTILS_C_

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

#endif /* UTILS_C_ */
