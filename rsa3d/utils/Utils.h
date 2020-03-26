/*
 * Utils.c
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */


#ifndef UTILS_C_
#define UTILS_C_

#include "../geometry/Vector.h"

template <std::size_t N>
using Orientation = std::array<double, N>;

using RSAVector = Vector<RSA_SPATIAL_DIMENSION>;
using RSAOrientation = Orientation<RSA_ANGULAR_DIMENSION>;

bool increment(int* in, int inlength, int max);
int position2i(const double* da, int dalength, double size, double dx, int n);
void i2position(double* da, int dalength, int index, double dx, int n);
void coordinates(int* result, const double* da, int dalength, double size, double dx, int n);
int neighbour2i(int* coordinates, int* neighbour, int clength, int offset, int n);

template<std::size_t SIZE>
void calculateGrid(std::array<int, SIZE> &grid, size_t index, size_t lenght){
	for(unsigned char i = 0; i<SIZE; i++){
		grid[i] = index % lenght;
		index /= lenght;
	}
}

// trim from start
std::string &ltrim(std::string &s);
// trim from end
std::string &rtrim(std::string &s);
// trim from both ends
std::string &trim(std::string &s);
// replaces all occurences of search in source by replace
std::string replaceAll(std::string source, const std::string& search, const std::string& replace);
bool endsWith(const std::string& str, const std::string& suffix);
bool startsWith(const std::string& str, const std::string& prefix);
int lastIndexOf(const std::string &s, char target);

void die(const std::string & reason);

double getAngleToOrigin(const Vector<2> &point);
void rotate2D(double* point, double alpha);

#endif /* UTILS_C_ */
