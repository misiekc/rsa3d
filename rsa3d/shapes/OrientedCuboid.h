/*
 * OrientedCuboid.h
 *
 *  Created on: 12.07.2017
 *      Author: ciesla
 */

#ifndef SHAPES_ORIENTEDCUBOID_H_
#define SHAPES_ORIENTEDCUBOID_H_


template <unsigned short DIMENSION>
class OrientedCuboid : public Shape<DIMENSION, 0>{
private:
//	static bool do2Drotation;
    static double           neighbourListCellSize;
    static double           voxelSize;

    static double 			size[DIMENSION];

public:
	OrientedCuboid();
	virtual ~OrientedCuboid();

	static void initClass(const std::string &args);
	static Shape<DIMENSION, 0> * create(RND *rnd);

    int overlap(BoundaryConditions *bc, Shape<DIMENSION, 0> *s);
    int pointInside(BoundaryConditions *bc, double* da);
	int pointInside(BoundaryConditions *bc, double* position, double *orientation, double orientationRange);
    double getNeighbourListCellSize();
    double getVoxelSize();
    double getVolume();

    std::string toPovray() const;
};

#include "OrientedCuboid.tpp"

#endif /* SHAPES_ORIENTEDCUBOID_H_ */
