//
// Created by Piotr Kubala on 26/03/2020.
//

#include <sstream>
#include <cmath>

#include "RegularDiskopolygon.h"
#include "spherocylinder/SpheroCylinder2D.h"
#include "../../FreeBC.h"

std::size_t RegularDiskopolygon::nSides{};
double RegularDiskopolygon::sideLength{};
double RegularDiskopolygon::radius{};
double RegularDiskopolygon::height{};
double RegularDiskopolygon::halfDiagonal{};

void RegularDiskopolygon::initClass(const std::string &attr) {
    std::istringstream attrStream(attr);
    attrStream >> nSides >> sideLength >> radius;
    ValidateMsg(attrStream, "Malformed attributes. Expected: [nubmer of sites] [side length] [disk radius]");
    Validate(nSides >= 3);
    Validate(sideLength >= 0);
    Validate(radius > 0);

    height = sideLength / 2 / std::tan(M_PI / nSides);
    halfDiagonal = sideLength / 2 / std::sin(M_PI / nSides);
    normalizeVolume();

    std::ostringstream spherocylinderAttr;
    spherocylinderAttr << sideLength << " " << radius;
    SpheroCylinder2D::initClass(spherocylinderAttr.str());

    ShapeStaticInfo<2, 1> shapeInfo;

    shapeInfo.setCircumsphereRadius(halfDiagonal + radius);
    shapeInfo.setInsphereRadius(height + radius);
    shapeInfo.setAngularVoxelSize(2*M_PI / nSides);
    shapeInfo.setSupportsSaturation(true);
    shapeInfo.setDefaultCreateShapeImpl<RegularDiskopolygon>();

    Shape::setShapeStaticInfo(shapeInfo);
}

void RegularDiskopolygon::normalizeVolume() {
    double volume = 0.25 * nSides * std::pow(sideLength, 2) / tan(M_PI / nSides)   // Polygon inside
                    + nSides * radius * sideLength                                  // Pushed edges
                    + M_PI * std::pow(radius, 2);                                   // Corners

    sideLength /= std::sqrt(volume);
    radius /= std::sqrt(volume);
    height /= std::sqrt(volume);
    halfDiagonal /= std::sqrt(volume);
}

bool RegularDiskopolygon::overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const {
    switch (this->overlapEarlyRejection(bc, s)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }

    RegularDiskopolygon other = dynamic_cast<const RegularDiskopolygon &>(*s);
    this->applyBC(bc, &other);

    FreeBC<2> freeBC;
    for (std::size_t i{}; i < nSides; i++) {
        SpheroCylinder2D sc1 = this->getSpherocylinder(i);
        for (std::size_t j{}; j < nSides; j++) {
            SpheroCylinder2D sc2 = other.getSpherocylinder(j);
            if (sc1.overlap(&freeBC, &sc2))
                return true;
        }
    }

    return false;
}

bool RegularDiskopolygon::voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition,
                                      const Orientation<1> &orientation, double spatialSize, double angularSize) const
{
    switch(this->voxelInsideEarlyRejection(bc, voxelPosition, orientation, spatialSize, angularSize)) {
        case EarlyRejectionResult::TRUE:      return true;
        case EarlyRejectionResult::FALSE:     return false;
        case EarlyRejectionResult::UNKNOWN:   break;
    }

    Vector<2> defaultSpherocylinderOffset{{0, height}};

    for (std::size_t i{}; i < nSides; i++) {
        double virtualScAngleFrom = orientation[0] + static_cast<double>(i) * 2 * M_PI / nSides;
        double virtualScAngleTo = virtualScAngleFrom + angularSize;
        AnisotropicShape2D::normalizeAngleRange(0, &virtualScAngleFrom, &virtualScAngleTo, 2*M_PI);

        RectangularBounding rectangularBounding = RectangularBoundingBuilder::buildForArch(defaultSpherocylinderOffset,
                virtualScAngleFrom, virtualScAngleTo);

        Vector<2> voxelTranslation = bc->getTranslation(this->getPosition(), voxelPosition);
        rectangularBounding.expand(spatialSize);
        rectangularBounding.translate(voxelPosition + voxelTranslation);

        for (std::size_t j{}; j < nSides; j++) {
            FreeBC<2> freeBC;

            auto sc = this->getSpherocylinder(j);
            if (sc.pointInside(&freeBC, rectangularBounding.getTopLeft(), virtualScAngleFrom, virtualScAngleTo) &&
                sc.pointInside(&freeBC, rectangularBounding.getTopRight(), virtualScAngleFrom, virtualScAngleTo) &&
                sc.pointInside(&freeBC, rectangularBounding.getBottomLeft(), virtualScAngleFrom, virtualScAngleTo) &&
                sc.pointInside(&freeBC, rectangularBounding.getBottomRight(), virtualScAngleFrom, virtualScAngleTo))
            {
                return true;
            }
        }
    }

    return false;
}

Shape<2, 1> *RegularDiskopolygon::clone() const {
    return new RegularDiskopolygon(*this);
}

SpheroCylinder2D RegularDiskopolygon::getSpherocylinder(std::size_t index) const {
    SpheroCylinder2D sc;

    double scAngle = this->getAngle() + static_cast<double>(index) * 2 * M_PI / nSides;
    Vector<2> defaultSpherocylinderOffset{{0, height}};
    sc.rotate({{scAngle}});
    sc.translate(this->getPosition() + Matrix<2, 2>::rotation(scAngle) * defaultSpherocylinderOffset);

    return sc;
}

std::string RegularDiskopolygon::toWolfram() const {
    std::ostringstream out;

    out << "{";
    for (std::size_t i{}; i < nSides; i++)
        out << this->getSpherocylinder(i).toWolfram() << ", ";
    // For this->getAngle() == 0, we want top side to be horizontal
    double mathematicaAngle = this->getAngle() + M_PI/2 + M_PI/nSides;
    out << "RegularPolygon[" << this->getPosition() << ", {" << halfDiagonal << ", " << mathematicaAngle << "}, ";
    out << nSides << "]}";

    return out.str();
}

void RegularDiskopolygon::setOrientation(const Orientation<1> &orientation) {
    double angle = AnisotropicShape2D::normalizeAngle(orientation[0], 2*M_PI/nSides);
    Shape::setOrientation({{angle}});
}

void RectangularBoundingBuilder::createBounding(RectangularBounding &bounding, const Vector<2> &zeroAngleVector,
                                                double angleTo, double quarterAngle)
{
    if (angleTo < quarterAngle) {
        bounding.addPoint(Matrix<2, 2>::rotation(angleTo) * zeroAngleVector);
    } else {
        bounding.addPoint(Matrix<2, 2>::rotation(quarterAngle) * zeroAngleVector);
        createBounding(bounding, zeroAngleVector, angleTo, quarterAngle + M_PI / 2);
    }
}

void RectangularBounding::addPoint(const Vector<2> &p) {
    if (p[0] < minPoint[0])
        minPoint[0] = p[0];
    if (p[1] < minPoint[1])
        minPoint[1] = p[1];
    if (p[0] > maxPoint[0])
        maxPoint[0] = p[0];
    if (p[1] > maxPoint[1])
        maxPoint[1] = p[1];
}

void RectangularBounding::translate(const Vector<2> &translation) {
    this->minPoint += translation;
    this->maxPoint += translation;
}

void RectangularBounding::expand(double expansion) {
    this->maxPoint += {{expansion, expansion}};
}

RectangularBounding RectangularBoundingBuilder::buildForArch(const Vector<2> &zeroAngleVector, double angleFrom,
                                                             double angleTo)
{
    Expects(angleFrom >= 0);
    Expects(angleTo >= angleFrom);

    RectangularBounding rectangularBounding;
    rectangularBounding.addPoint(Matrix<2, 2>::rotation(angleFrom) * zeroAngleVector);

    double quarterAngle = M_PI/2;
    while (quarterAngle < angleFrom)
        quarterAngle += M_PI/2;

    createBounding(rectangularBounding, zeroAngleVector, angleTo, quarterAngle);
    return rectangularBounding;
}
