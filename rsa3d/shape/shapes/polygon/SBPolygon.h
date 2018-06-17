/*
 * SBPolygon.h
 *
 *  Created on: 04.06.2018
 *      Author: ciesla
 */

#ifndef SHAPES_POLYGONS_SBPOLYGON_H_
#define SHAPES_POLYGONS_SBPOLYGON_H_

#include "Polygon.h"

class SBPolygon: public Polygon {
public:
	SBPolygon();
	virtual ~SBPolygon();

	static void initClass(const std::string &args);

};

#endif /* SHAPES_POLYGONS_SBPOLYGON_H_ */
