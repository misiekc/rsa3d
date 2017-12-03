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
#include "../Vector.h"


class OverlapStrategy;

class Cuboid : public Shape<3>
{
private:
    static double           volume;
    static double           *auxDoubleArray;        // Auxiliary double array of dimension size
    Matrix<3, 3>            orientation;
    static double           minDimension;
    static OverlapStrategy  *strategy;
    static Vector<3>        relativeVertices[8];

protected:
    static unsigned short 	staticDimension;
    static double           *size;
    
    static double           neighbourListCellSize;
    static double           voxelSize;

public:

    // Vertex recognition helper. P states positive, N - negative. First position
    // corresponds to positive/negative X, second for Y, etc.
    enum VERTEX {
        PPP = 0,
        NPP,
        PNP,
        PPN,
        PNN,
        NPN,
        NNP,
        NNN,
        NUM_OF
    };

    // Coord recognition helper
    enum COORD {
        X = 0,
        Y,
        Z
    };

    // Implicit copy ctor and copy assignment operator - trivial destructor
    explicit Cuboid(const Matrix<3, 3> & rotation);
    ~Cuboid() override = default;

	static void initClass(const std::string &args);
	static Shape<3> * create2D(RND *rnd);
	static Shape<3> * create3D(RND *rnd);
	static void setOverlapStrategy(OverlapStrategy *strategy);
	static OverlapStrategy * getOverlapStrategy();
    static const Vector<3> getRelativeVertex(std::size_t index);

    double getNeighbourListCellSize() override;
    double getVoxelSize() override;
    int overlap(BoundaryConditions *bc, Shape<3> *s) override;
    double getVolume() override;
    int pointInside(BoundaryConditions *bc, double* da) override;

    bool pointInsideCuboid(const Vector<3> &vertex);
    void obtainVertices(Vector<3> (&vertices)[8], const Vector<3> &translation);
    
    Matrix<3, 3> getOrientation() const;
    static double * getSize(double * arr);
    std::string toPovray() const override;
    std::string toWolfram() const;
	void store(std::ostream &f) const override;
	void restore(std::istream &f) override;
};


#endif      // _CUBOID_H
