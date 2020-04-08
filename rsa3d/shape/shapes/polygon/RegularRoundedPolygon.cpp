//
// Created by pkua on 06.04.2020.
//

#include "RegularRoundedPolygon.h"
#include "../RegularDiskopolygon.h"

void RegularRoundedPolygon::initClass(const std::string &attr) {
    RegularDiskopolygonAttributes attributes(attr);

    double nSides = attributes.getNSides();
    double radius = attributes.getRadius();
//    double height = attributes.getHeight();
    double halfDiagonal = attributes.getHalfDiagonal();

    std::ostringstream roundedPolygonAttributes;
    roundedPolygonAttributes.precision(std::numeric_limits< double >::max_digits10);

    roundedPolygonAttributes << radius << " " << nSides << " rt ";

    for (std::size_t i{}; i < nSides; i++)
//        roundedPolygonAttributes << " " << height << " " << (2*M_PI*i/nSides);
        roundedPolygonAttributes << " " << halfDiagonal << " " << (2*M_PI*i/nSides);

    roundedPolygonAttributes << " " << nSides;
    for (std::size_t i{}; i < nSides; i++)
        roundedPolygonAttributes << " " << i;

    RoundedPolygon::initClass(roundedPolygonAttributes.str());

    ShapeStaticInfo<2, 1> shapeInfo = Shape::getShapeStaticInfo();
    shapeInfo.setAngularVoxelSize(2*M_PI/nSides);
    Shape::setShapeStaticInfo(shapeInfo);

}
