//----------------------------------------------------------------------------
// Class representing cuboid, derived after Shape class.
//----------------------------------------------------------------------------
// (C)PKua 2017
//----------------------------------------------------------------------------


#ifndef _CUBOID_H
#define _CUBOID_H

#define CUBOID_DEBUG        // Used for debbuging output

#include "../Shape.h"
#include "../Matrix.h"

class Cuboid : public Shape
{
private:
    static double           *size;
    static double           volume;
    static unsigned char    staticDimension;
    Matrix                  rotation;
    
    static double           neighbourListCellSize;
    static double           voxelSize;

public:
    Cuboid(unsigned char dimension, const Matrix & rotation);
    ~Cuboid();

	static void initClass(const std::string &args);
	static Shape * create(RND *rnd);

    double getNeighbourListCellSize();
    double getVoxelSize();
    int overlap(BoundaryConditions *bc, Shape *s);
    double getVolume();
    int pointInside(BoundaryConditions *bc, double* da);
};

#endif      // _CUBOID_H
