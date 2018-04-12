/*
 * Voxel.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef VOXEL_H_
#define VOXEL_H_

#include <algorithm>

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
class Voxel{

private:
	double position[SPATIAL_DIMENSION];
	std::array<double, ANGULAR_DIMENSION> orientation;
	short missCounter;

public:
	/**
	 * index of the voxel in array used by VoxelList - used for faster removing from the array
	 */
	int index;

	/**
	 * number of particle in the packing where voxel was checked if it is inside an exclusion zone - used for faster checking in the future
	 */
	int lastAnalyzed;

	/**
	 * depth level at which the voxel was checked if it is inside an exclusion zone
	 */
	unsigned short depth;


	/**
	 * creates an empty voxel. It position and orientation is not initialized and other attributes are set to 0
	 */
	Voxel();

	/**
	 * creates a voxel of a given position and orientation
	 */
	Voxel(double* pos, double *angle);

	virtual ~Voxel() = default;

	void miss();

	int getMissCounter();

	void resetMissCounter();

	bool isInside(double *pos, double size);

	bool isInside(double *pos, double size, const double* angle, double asize);

	double *getPosition();

	double *getOrientation();

	std::string toPovray(double ssize);
	std::string toWolfram(double ssize, double asize);
	std::string toString();

	void store(std::ostream &f) const;
	void restore(std::istream &f);



};

#include "Voxel.tpp"

#endif /* VOXEL_H_ */
