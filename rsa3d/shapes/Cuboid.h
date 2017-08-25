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
public:
    // Overlap stategy
    enum OverlapStrategy {
        MINE,
        TRI_TRI,
        SAT
    };

private:
    static double           volume;
    static double           *auxDoubleArray;        // Auxiliary double array of dimension size
    static double           *auxDoubleArray2;       // Second aux array
    Matrix<3, 3>            orientation;
    static double           minDimension;
    
    static OverlapStrategy  strategy;
    
    bool        checkPoint(const Vector<3> & vertex);
    bool        checkSegment(const Vector<3> & point1, const Vector<3> & point2);
    
    int overlapMine(BoundaryConditions *bc, Shape *s);
    int overlapTri(BoundaryConditions *bc, Shape *s);
    int overlapSAT(BoundaryConditions *bc, Shape *s);
    
    void obtainTris(Vector<3> (&arr)[12][3], const Vector<3> & translation);

protected:
    static double           *size;
    static unsigned char    staticDimension;
    
    static double           neighbourListCellSize;
    static double           voxelSize;

public:

    // Implicit copy ctor and copy assignment operator - trivial destructor
    Cuboid(const Matrix<3, 3> & rotation);
    ~Cuboid();

	static void initClass(const std::string &args);
	static Shape * create(RND *rnd);
	static void setOverlapStrategy(OverlapStrategy strategy);
	static OverlapStrategy getOverlapStrategy();

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
