//----------------------------------------------------------------------------
// Class representing cuboid, derived after Shape class.
//----------------------------------------------------------------------------
// (C)PKua 2017
//----------------------------------------------------------------------------


#ifndef _CUBOID_H
#define _CUBOID_H

//#define CUBOID_DEBUG        // Used for debbuging output

#include "../../../geometry/Matrix.h"
#include "../../../geometry/Vector.h"
#include "../../OverlapStrategyShape.h"
#include "../../../geometry/Geometry.h"
#include "../../ConvexShape.h"
#include "../../OrderCalculable.h"


class CuboidOverlapStrategy;

class Cuboid : public ConvexShape<3, 0>, public OverlapStrategyShape<3, 0>, public OrderCalculable
{
private:
	static double           size[3];
	static double           insphereRadius;
	static double           circumsphereRadius;
	
	static Vector<3>        relativeVertices[8];
	
	static CuboidOverlapStrategy  *defaultStrategy;
	static CuboidOverlapStrategy  *strategy;

    Matrix<3, 3>            orientation;

	static void calculateRelativeVerties();
    static Shape<3, 0> * create2D(RND *rnd);
    static Shape<3, 0> * create3D(RND *rnd);
    bool liesInCylinderOnEdge(const Vector<3> &absPointPos, size_t coord1, size_t coord2) const;

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
	static double * getSize(double * arr);
	static const Vector<3> getRelativeVertex(std::size_t index);

	static void restoreDefaultStrategy();
	static void setOverlapStrategy(CuboidOverlapStrategy *strategy);
	static CuboidOverlapStrategy * getOverlapStrategy();

	bool overlap(BoundaryConditions<3> *bc, const Shape<3,0> *s) const override;
	double getVolume() const override;
	bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     double orientationRange) const override;
	std::string toPovray() const override;
	std::string toWolfram() const override;
	void store(std::ostream &f) const override;
	Shape<3, 0> *clone() const override;

	void restore(std::istream &f) override;

	std::vector<std::string> getSupportedStrategies() const override;
	OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;

	std::vector<double> calculateOrder(const OrderCalculable *other) const override;

	bool pointInsideCuboid(const Vector<3> &vertex) const;
	void obtainVertices(Vector<3> (&vertices)[8], const Vector<3> &translation) const;
	PolyhedronTriangulation obtainTris() const;
	Matrix<3, 3> getOrientation() const;

    static void normalizeVolume();
};


#endif      // _CUBOID_H
