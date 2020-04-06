//
// Created by pkua on 06.04.2020.
//

#include "RegularRoundedPolygon.h"
#include "../RegularDiskopolygon.h"

void RegularRoundedPolygon::initClass(const std::string &attr) {
    RegularDiskopolygonAttributes attributes(attr);

    double nSides = attributes.getNSides();
    double radius = attributes.getRadius();
    double height = attributes.getHeight();

    std::ostringstream roundedPolygonAttributes;
    roundedPolygonAttributes << radius << " " << nSides << " rt ";

    for (std::size_t i{}; i < nSides; i++)
        roundedPolygonAttributes << " " << height << " " << (2*M_PI*i/nSides);

    roundedPolygonAttributes << " " << nSides;
    for (std::size_t i{}; i < nSides; i++)
        roundedPolygonAttributes << " " << i;

    RoundedPolygon::initClass(roundedPolygonAttributes.str());
}
