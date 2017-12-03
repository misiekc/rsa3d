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

#include <utility>


class OverlapStrategy;

class Cuboid : public Shape<3>
{
public:


private:
    static double           volume;
    static double           *auxDoubleArray;        // Auxiliary double array of dimension size
    static double           *auxDoubleArray2;       // Second aux array
    Matrix<3, 3>            orientation;
    static double           minDimension;
    static OverlapStrategy  *strategy;

    int overlapTri(BoundaryConditions *bc, Shape<3> *s);
    int overlapSAT(BoundaryConditions *bc, Shape<3> *s);
    
    void obtainTris(Vector<3> (&arr)[12][3], const Vector<3> & translation);

protected:
    static unsigned short 	staticDimension;
    static double           *size;
    
    static double           neighbourListCellSize;
    static double           voxelSize;

public:

    // Implicit copy ctor and copy assignment operator - trivial destructor
    explicit Cuboid(const Matrix<3, 3> & rotation);
    ~Cuboid() override;

	static void initClass(const std::string &args);
	static Shape<3> * create2D(RND *rnd);
	static Shape<3> * create3D(RND *rnd);
	static void setOverlapStrategy(OverlapStrategy *strategy);
	static OverlapStrategy * getOverlapStrategy();

    double getNeighbourListCellSize() override;
    double getVoxelSize() override;
    int overlap(BoundaryConditions *bc, Shape<3> *s) override;
    double getVolume() override;
    int pointInside(BoundaryConditions *bc, double* da) override;

    bool checkPoint(const Vector<3> & vertex);
    
    Matrix<3, 3> getOrientation() const;
    static double * getSize(double * arr);
    std::string toPovray() const override;
    std::string toWolfram() const;
	void store(std::ostream &f) const override;
	void restore(std::istream &f) override;
};


// Overlap stategies
class OverlapStrategy {
public:
    virtual bool overlap(Cuboid * cube1, Cuboid * cube2, BoundaryConditions *bc) = 0;
    virtual std::string getName() = 0;
};

class MineOverlap : public OverlapStrategy {
private:
    bool checkSegment(Cuboid *cube, const Vector<3> & point1, const Vector<3> & point2);

public:
    bool overlap(Cuboid *cube1, Cuboid *cube2, BoundaryConditions *bc) override;
    std::string getName() override;
};

class SATOverlap : public OverlapStrategy {
private:
    typedef std::pair<double, double>   interval;

    bool checkSeparatingAxis(const Vector<3> & _axis, Vector<3> * _vert1, Vector<3> * _vert2) const;
    interval getProjection(const Vector<3> & _axis, Vector<3> * _vert) const;

public:
    bool overlap(Cuboid *cube1, Cuboid *cube2, BoundaryConditions *bc) override;
    std::string getName() override;
};

#endif      // _CUBOID_H
