/*
 * Voxel.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef VOXEL_H_
#define VOXEL_H_

#include <algorithm>
#include <array>
#include "Vector.h"

class Voxel{

private:
	using RSAVector = Vector<RSA_SPATIAL_DIMENSION>;

	RSAVector position;
	std::array<double, RSA_ANGULAR_DIMENSION> orientation;

public:

	/**
	 * @brief number of particle in the packing where voxel was checked if it is inside an exclusion zone - used for faster checking in the future
	 */
	int lastAnalyzed;

	/**
	 * @brief depth level at which the voxel was checked if it is inside an exclusion zone
	 */
	unsigned short depth;


	/**
	 * @brief creates an empty voxel. It position and orientation is not initialized and other attributes are set to 0
	 */
	Voxel();

	/**
	 * @brief creates a voxel of a given position and orientation
	 */
	Voxel(const RSAVector &pos, const std::array<double, RSA_ANGULAR_DIMENSION> &angle);

	virtual ~Voxel() = default;

	bool isInside(const RSAVector &pos, double size);

	bool isInside(const RSAVector &pos, double size, const std::array<double, RSA_ANGULAR_DIMENSION> &angle,
                  double asize);

	const RSAVector &getPosition();

	std::array<double, RSA_ANGULAR_DIMENSION> getOrientation();

	std::string toPovray(double ssize);
	std::string toWolfram(double ssize, double asize);
	std::string toString();

	void store(std::ostream &f) const;
	void restore(std::istream &f);
};


#endif /* VOXEL_H_ */
