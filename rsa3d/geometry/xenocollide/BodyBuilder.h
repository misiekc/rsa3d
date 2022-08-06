/*
 * BodyBuilder.h
 *
 *  Created on: Aug 5, 2022
 *      Author: ciesla
 */

#ifndef GEOMETRY_XENOCOLLIDE_BODYBUILDER_H_
#define GEOMETRY_XENOCOLLIDE_BODYBUILDER_H_

#include <list>
#include "CollideGeometry.h"
#include "Quat.h"
#include "MapPtr.h"
#include "../Vector.h"

class BodyBuilder
{

public:

	BodyBuilder();
	~BodyBuilder();

	// shapes
	void axis(double x);
	void box(double x, double y, double z);
	void disc(double x);
	void disc(double x, double y);
	void football(double l, double w);
	void point(double x, double y, double z);
	void rect(double x, double y);
	void saucer(double r, double t);
	void segment(double l);
	void sphere(double x);
	void sphere(double x, double y, double z);

	// shapes transformations
	void move(double x, double y, double z);
	void rot(double x, double y, double z);

	// shape combination
	void diff();
	void sweep();
	void wrap();

	// stack operations
	void dup(size_t n);
	void pop();
	void swap();
	void clear();

	void ProcessCommand(std::string& cmd);
	MapPtr<CollideGeometry> getCollideGeometry();
	double getMaxRadius();

private:
	struct XCShape{
		XCShape() : geom(NULL) { q = Quat(0, 0, 0, 1); x = Vector<3>({0, 0, 0}); }
		XCShape(CollideGeometry* _geom, const Quat& _q, const Vector<3>& _x) : geom(_geom), q(_q), x(_x) {}
		XCShape(CollideGeometry* _geom) : geom(_geom) { q = Quat(0, 0, 0, 1); x = Vector<3>({0, 0, 0}); }
		XCShape(const XCShape& s) : geom(s.geom), q(s.q), x(s.x) {}
		XCShape& operator = (const XCShape& s){ geom = s.geom; q = s.q; x = s.x; return *this;}
		MapPtr<CollideGeometry>	geom;
		Quat				q;
		Vector<3>			x;
	};
	std::list< MapPtr<XCShape> > mShapeStack;
	std::string mStatusLine;
	void Status(const char* str);


};


#endif /* GEOMETRY_XENOCOLLIDE_BODYBUILDER_H_ */
