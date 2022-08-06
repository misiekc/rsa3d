/*
XenoCollide Collision Detection and Physics Library
Copyright (c) 2007-2014 Gary Snethen http://xenocollide.com

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising
from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must
not claim that you wrote the original software. If you use this
software in a product, an acknowledgment in the product documentation
would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must
not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
/*
 * Adapted by Michał Cieśla
 */

#include "Collide.h"
#include "CollideGeometry.h"
#include "../Vector.h"
#include "../../utils/Assertions.h"

//////////////////////////////////////////////////////////////////////////////

inline Vector<3> Collide::TransformSupportVert( CollideGeometry& p, const Quat& q, const Vector<3>& t, const Vector<3>& n )
{
	Vector<3> localNormal = (~q).Rotate(n);
	Vector<3> localSupport = p.GetSupportPoint(localNormal);
	Vector<3> worldSupport = q.Rotate(localSupport) + t;
	return worldSupport;
}

//////////////////////////////////////////////////////////////////////////////
// Collide

// Global Test Data
double gAvgSupportCount = 0.0f;
int gMaxPhase1 = 0;
int gMaxPhase2 = 0;
int gMinPhase1 = 1000;
int gMinPhase2 = 1000;
double gAvgPhase1 = 0;
double gAvgPhase2 = 0;
double gCountPhase1 = 0;
double gCountPhase2 = 0;

//////////////////////////////////////////////////////////////////////////////

inline void Collide::Swap(Vector<3>& a, Vector<3>& b)
{
	Vector<3> tmp = a;
	a = b;
	b = tmp;
}

///
// bool Intersect(...)
//
// Returns true if two CollideGeometry objects intersect.
//
bool Collide::Intersect(CollideGeometry& p1, const Quat& q1, const Vector<3>& t1, CollideGeometry& p2, const Quat& q2, const Vector<3>& t2, double boundaryTolerance)
{
	// v0 = center of Minkowski difference
	Vector<3> v0 = (q2.Rotate(p2.GetCenter()) + t2) - (q1.Rotate(p1.GetCenter()) + t1);
	if (v0.isZero()) return true;	// v0 and origin overlap ==> hit

	// v1 = support in direction of origin
	Vector<3> n = -v0;
	Vector<3> v1 = TransformSupportVert(p2, q2, t2, n) - TransformSupportVert(p1, q1, t1, -n);
	if (v1 * n <= 0) return false;	// origin outside v1 support plane ==> miss

	// v2 = support perpendicular to plane containing origin, v0 and v1
	n = v1 ^ v0;
	if (n.isZero()) return true;	// v0, v1 and origin colinear (and origin inside v1 support plane) == > hit
	Vector<3> v2 = TransformSupportVert(p2, q2, t2, n) - TransformSupportVert(p1, q1, t1, -n);
	if (v2 * n <= 0) return false;	// origin outside v2 support plane ==> miss

	// v3 = support perpendicular to plane containing v0, v1 and v2
	n = (v1 - v0) ^ (v2 - v0);

	// If the origin is on the - side of the plane, reverse the direction of the plane
	if (n * v0 > 0)
	{
		Collide::Swap(v1, v2);
		n = -n;
	}

	///
	// Phase One: Find a valid portal

	while (1)
	{
		// Obtain the next support point
		Vector<3> v3 = TransformSupportVert(p2, q2, t2, n) - TransformSupportVert(p1, q1, t1, -n);
		if (v3 * n <= 0) return false;	// origin outside v3 support plane ==> miss

		// If origin is outside (v1,v0,v3), then portal is invalid -- eliminate v2 and find new support outside face
		if ( (v1 ^ v3) * v0 < 0)
		{
			v2 = v3;
			n = (v1 - v0) ^ (v3 - v0);
			continue;
		}

		// If origin is outside (v3,v0,v2), then portal is invalid -- eliminate v1 and find new support outside face
		if ( (v3 ^ v2) * v0 < 0)
		{
			v1 = v3;
			n = (v3 - v0) ^ (v2 - v0);
			continue;
		}

		///
		// Phase Two: Refine the portal

		while (1)
		{
			// Compute outward facing normal of the portal
			n = (v2 - v1) ^ (v3 - v1);

			// If the origin is inside the portal, we have a hit
			if (n * v1 >= 0) return true;

			// Find the support point in the direction of the portal's normal
			Vector<3> v4 = TransformSupportVert(p2, q2, t2, n) - TransformSupportVert(p1, q1, t1, -n);

			// If the origin is outside the support plane or the boundary is thin enough, we have a miss
			n = n.normalized();
			if ( -(v4 * n) >= 0 || (v4 - v3) * n <= boundaryTolerance ) return false;

#if 1
			// Test origin against the three planes that separate the new portal candidates: (v1,v4,v0) (v2,v4,v0) (v3,v4,v0)
			// Note:  We're taking advantage of the triple product identities here as an optimization
			//        (v1 % v4) * v0 == v1 * (v4 % v0)    > 0 if origin inside (v1, v4, v0)
			//        (v2 % v4) * v0 == v2 * (v4 % v0)    > 0 if origin inside (v2, v4, v0)
			//        (v3 % v4) * v0 == v3 * (v4 % v0)    > 0 if origin inside (v3, v4, v0)
			Vector<3> cross = v4 ^ v0;
			if (v1 * cross > 0)
			{
				if (v2 * cross > 0) v1 = v4;	// Inside v1 & inside v2 ==> eliminate v1
				else v3 = v4;					// Inside v1 & outside v2 ==> eliminate v3
			}
			else
			{
				if (v3 * cross > 0) v2 = v4;	// Outside v1 & inside v3 ==> eliminate v2
				else v1 = v4;					// Outside v1 & outside v3 ==> eliminate v1
			}
#else
			// Test origin against the three planes that separate the new portal candidates: (v1,v4,v0) (v2,v4,v0) (v3,v4,v0)
			// Note: non-optimized version
			if (v1 ^ v4 * v0 > 0)
			{
				if (v2 ^ v4 * v0 > 0) v1 = v4;	// Inside v1 & inside v2 ==> eliminate v1
				else v3 = v4;					// Inside v1 & outside v2 ==> eliminate v3
			}
			else
			{
				if (v3 ^ v4 * v0 > 0) v2 = v4;	// Outside v1 & inside v3 ==> eliminate v2
				else v1 = v4;					// Outside v1 & outside v3 ==> eliminate v1
			}
#endif
		}
	}
}
