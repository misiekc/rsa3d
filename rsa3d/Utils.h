/*
 * Utils.c
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */


#ifndef UTILS_C_
#define UTILS_C_

bool increment(int* in, int inlength, int max);
int position2i(double* da, int dalength, double size, double dx, int n);
void coordinates(int* result, double* da, int dalength, double size, double dx, int n);
int neighbour2i(int* coordinates, int* neighbour, int clength, int offset, int n);

#endif /* UTILS_C_ */
