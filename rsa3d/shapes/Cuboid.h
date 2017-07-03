//----------------------------------------------------------------------------
// Class representing cuboid, derived after Shape class.
//----------------------------------------------------------------------------
// (C)PKua 2017
//----------------------------------------------------------------------------


#ifndef _CUBOID_H
#define _CUBOID_H

//#define CUBOID_DEBUG        // Used for debbuging output

#include "../Shape.h"
#include "../Matrix.h"

class Cuboid : public Shape
{
private:
    static double           *size;
    static double           volume;
    static unsigned char    staticDimension;
    static double           *auxDoubleArray;        // Auxiliary double array of dimension size
    static double           *auxDoubleArray2;       // Second aux array
    Matrix                  orientation;
    static double           minDimension;
    
    static double           neighbourListCellSize;
    static double           voxelSize;
    
    bool        checkPoint(const Matrix & vertex);
    bool        checkSegment(const Matrix & point1, const Matrix & point2);

public:
    // Implicit copy ctor and copy assignment operator - trivial destructor
    Cuboid(const Matrix & rotation);
    ~Cuboid();

	static void initClass(const std::string &args);
	static Shape * create(RND *rnd);

    double getNeighbourListCellSize();
    double getVoxelSize();
    int overlap(BoundaryConditions *bc, Shape *s);
    double getVolume();
    int pointInside(BoundaryConditions *bc, double* da);
    
    Matrix getOrientation() const;
};

#endif      // _CUBOID_H
