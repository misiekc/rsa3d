//
// Created by Michal Ciesla on 10.01.2023.
//

#include "RegularPolygon.h"
#include "../RegularDiskopolygon.h"

void RegularPolygon::initClass(const std::string &attr) {
    RegularDiskopolygonAttributes attributes(attr);

    double nSides = attributes.getNSides();
    double halfDiagonal = attributes.getHalfDiagonal();

    std::ostringstream polygonAttributes;
    polygonAttributes.precision(std::numeric_limits< double >::max_digits10);

    polygonAttributes << nSides << " rt ";

    for (std::size_t i{}; i < nSides; i++)
        polygonAttributes << " " << halfDiagonal << " " << (2*M_PI*i/nSides);

    polygonAttributes << " " << nSides;
    for (std::size_t i{}; i < nSides; i++)
        polygonAttributes << " " << i;

    polygonAttributes << " starHelperSegments";

    Polygon::initClass(polygonAttributes.str());

    // Make angular voxel size smaller - regular polygon has an n-fold rotational symmetry
    ShapeStaticInfo<2, 1> shapeInfo = Shape::getShapeStaticInfo();
    shapeInfo.setAngularVoxelSize(2*M_PI/nSides);
    Shape::setShapeStaticInfo(shapeInfo);
}
