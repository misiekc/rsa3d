//
// Created by Kasperek, Wojciech on 10/03/2018.
//

#ifndef RSA3D_RECTANGLE_H
#define RSA3D_RECTANGLE_H

#include "../Shape.h"
#include "../RND.h"
#include "../BoundaryConditions.h"
#include "../AnisotropicShape2D.h"
#include "math.h"
#include "../Vector.h"

class Rectangle: public AnisotropicShape2D {
private:
    static double longer;
    static double shorter;

    // b is shorter
    double a, b, halfA, halfB;

    // vertices, where last = first
    double xs[5];
    double ys[5];

    double getXOnCircle(double cx, double cy, double r, double angle) const;

    double getYOnCircle(double cx, double cy, double r, double angle) const;

    int isInsideCircle(const Vector<2, double> &point, double cx, double cy, double r) const;

    int isInsideRect(const Vector<2, double> &point, double cx, double cy, double aLength, double bLength) const;

    int pointInsideInternal(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) const;

    Vector<2, double> rotatePoint(const Vector<2, double>& point, const Matrix<2, 2, double>& rotation, const Vector<2, double>& center) const;

    Vector<2> rotatePoint(const Vector<2>& point, const Matrix<2,2>& rotation) const;

    bool isInsideExcludingRectangle(double angleFrom, double angleTo, double x, double y) const;

    bool checkTangents(double x, double y, double cx, double cy, double r, double angle1, double angle2) const;

    void setAngle(double angle) override;

    bool isIntersection(double p0_x, double p0_y, double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y) const;
    static Shape<2, 1> * create(RND *rnd);

protected:
    void setPosition(const double *position) override;

public:
    static void initClass(const std::string &args);

    Rectangle();
    ~Rectangle() override = default;

    int overlap(BoundaryConditions *bc, Shape<2, 1> *s) const override;

    double getVolume() const override;

    int pointInside(BoundaryConditions *bc, double* da) const override;
    int pointInside(BoundaryConditions *bc, double* da, double angleFrom, double angleTo) const override;

    std::string toWolfram() const override;
    std::string toPovray() const override;
    std::string toString() const override;

    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
};


#endif //RSA3D_RECTANGLE_H
