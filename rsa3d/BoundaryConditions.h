/*
 * BoundaryConditions.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef BOUNDARYCONDITIONS_H_
#define BOUNDARYCONDITIONS_H_

#include "Vector.h"

template <unsigned short SPATIAL_DIMENSION>
class BoundaryConditions {
public:
    using vector = Vector<SPATIAL_DIMENSION>;

	virtual ~BoundaryConditions() = default;

	virtual double distance2(const vector &p1, const vector &p2) const = 0;

	/**
	 * @brief Returns translation that should be applied to @a p2 to move him to the "proximity" of @a p1
	 */
	virtual vector getTranslation(const vector& p1, const vector &p2) const = 0;

};

using RSABoundaryConditions = BoundaryConditions<RSA_SPATIAL_DIMENSION>;

#endif /* BOUNDARYCONDITIONS_H_ */
