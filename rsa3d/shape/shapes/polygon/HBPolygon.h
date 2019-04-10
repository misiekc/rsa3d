/*
 * HBPolygon.h
 *
 *  Created on: 04.06.2018
 *      Author: ciesla
 */

#ifndef SHAPES_POLYGONS_HBPOLYGON_H_
#define SHAPES_POLYGONS_HBPOLYGON_H_

#include "Polygon.h"

class HBPolygon : public Polygon{
public:
	~HBPolygon() override = default;

	static void initClass(const std::string &args);
};

#endif /* SHAPES_POLYGONS_HBPOLYGON_H_ */
