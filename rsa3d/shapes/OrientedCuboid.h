/*
 * OrientedCuboid.h
 *
 *  Created on: 12.07.2017
 *      Author: ciesla
 */

#ifndef SHAPES_ORIENTEDCUBOID_H_
#define SHAPES_ORIENTEDCUBOID_H_


template <unsigned short DIMENSION>
class OrientedCuboid : public Shape<DIMENSION>{
private:
//	static bool do2Drotation;
    static double           neighbourListCellSize;
    static double           voxelSize;

    static double 			size[DIMENSION];

public:
	OrientedCuboid();
	virtual ~OrientedCuboid();

	static void initClass(const std::string &args);
	static Shape<DIMENSION> * create(RND *rnd);

    int overlap(BoundaryConditions *bc, Shape<DIMENSION> *s);
    int pointInside(BoundaryConditions *bc, double* da);
    double getNeighbourListCellSize();
    double getVoxelSize();
    double getVolume();

    std::string toPovray() const;
};

#include "OrientedCuboid.tpp"

#endif /* SHAPES_ORIENTEDCUBOID_H_ */
