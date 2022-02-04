/*
 * RoundedPolygon.h
 *
 *  Created on: 3.11.2021
 *      Author: ciesla
 */

#ifndef SHAPES_POLYGONS_ROUNDEDRECTANGLE_H_
#define SHAPES_POLYGONS_ROUNDEDRECTANGLE_H_

#include "RoundedPolygon.h"

class RoundedRectangle : public RoundedPolygon{

private:
	std::vector<Vector<2>> getAxes();

public:
	static void initClass(const std::string &args);
    std::vector<double> calculateOrder(const OrderCalculable *other) const override;
};



#endif /* SHAPES_POLYGONS_ROUNDEDRECTANGLE_H_ */
