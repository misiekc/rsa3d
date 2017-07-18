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
    Matrix<3, 3>            orientation;
    static double           minDimension;
    
    static double           neighbourListCellSize;
    static double           voxelSize;
    
    bool        checkPoint(const Vector<3> & vertex);
    bool        checkSegment(const Vector<3> & point1, const Vector<3> & point2);

public:
    // Implicit copy ctor and copy assignment operator - trivial destructor
    Cuboid(const Matrix<3, 3> & rotation);
    ~Cuboid();

	static void initClass(const std::string &args);
	static Shape * create(RND *rnd);

    double getNeighbourListCellSize();
    double getVoxelSize();
    int overlap(BoundaryConditions *bc, Shape *s);
    double getVolume();
    int pointInside(BoundaryConditions *bc, double* da);
    
    Matrix<3, 3> getOrientation() const;
    static double * getSize(double * arr);
    std::string toPovray() const;
    std::string toWolfram() const;
	void store(std::ostream &f) const;
	void restore(std::istream &f);
};

#endif      // _CUBOID_H
