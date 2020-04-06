//
// Created by Piotr Kubala on 26/03/2020.
//

#ifndef RSA3D_REGULARDISKOPOLYGON_H
#define RSA3D_REGULARDISKOPOLYGON_H

#include "../Shape.h"
#include "spherocylinder/SpheroCylinder2D.h"

class RectangularBounding {
private:
    Vector<2> minPoint{{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()}};
    Vector<2> maxPoint{{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()}};

public:
    Vector<2> getBottomLeft() const { return minPoint; }
    Vector<2> getBottomRight() const { return {{maxPoint[0], minPoint[1]}}; }
    Vector<2> getTopRight() const { return maxPoint; }
    Vector<2> getTopLeft() const { return {{minPoint[0], maxPoint[1]}}; }

    void addPoint(const Vector<2> &p);
    void translate(const Vector<2> &translation);
    void expand(double expansion);
};

class RectangularBoundingBuilder {
private:
    static void recurseQuarters(RectangularBounding &bounding, const Vector<2> &zeroAngleVector, double angleTo,
                                double quarterAngle);

public:
    static RectangularBounding forArch(const Vector<2> &zeroAngleVector, double angleFrom, double angleTo);
};

class RegularDiskopolygonAttributes {
private:
    std::size_t nSides{};
    double sideLength{};
    double radius{};
    double height{};
    double halfDiagonal{};

    void normalizeVolume();

public:
    RegularDiskopolygonAttributes() = default;

    /**
     * @brief Calculates attributes from format: [number of sides] [side length] [radius]
     */
    explicit RegularDiskopolygonAttributes(const std::string &attr);

    size_t getNSides() const { return nSides;}
    double getSideLength() const { return sideLength; }
    double getRadius() const { return radius; }
    double getHeight() const { return height; }
    double getHalfDiagonal() const { return halfDiagonal; }
};

class RegularDiskopolygon : public Shape<2, 1> {
private:
    static RegularDiskopolygonAttributes attributes;

    SpheroCylinder2D getSpherocylinder(std::size_t index) const;

protected:
    void setOrientation(const Orientation<1> &orientation) override;

public:
    static void initClass(const std::string &attr);

    static std::size_t getNSides() { return attributes.getNSides(); }
    static double getSideLength() { return attributes.getSideLength(); }
    static double getRadius() { return attributes.getRadius(); }
    static double getHeight() { return attributes.getHeight(); }
    static double getHalfDiagonal() { return attributes.getHalfDiagonal(); }

    bool overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const override;
    bool voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition, const Orientation<1> &orientation,
                     double spatialSize, double angularSize) const override;
    Shape<2, 1> *clone() const override;
    std::string toWolfram() const override;
    double getVolume(unsigned short dim) const override { return 1; }

    double getAngle() const { return this->getOrientation()[0]; };
};


#endif //RSA3D_REGULARDISKOPOLYGON_H
