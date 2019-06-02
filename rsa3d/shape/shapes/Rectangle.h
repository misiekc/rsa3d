//
// Created by Kasperek, Wojciech on 10/03/2018.
//

#ifndef RSA3D_RECTANGLE_H
#define RSA3D_RECTANGLE_H

#include "../Shape.h"
#include "../../RND.h"
#include "../../BoundaryConditions.h"
#include "../AnisotropicShape2D.h"
#include "../OrderCalculable.h"
#include "../../geometry/Vector.h"
#include <cmath>

class Rectangle: public AnisotropicShape2D, public OrderCalculable {
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
    bool pointInsideInternal(BoundaryConditions<2> *bc, const Vector<2> &da, double angleFrom, double angleTo) const;
    Vector<2, double> rotatePoint(const Vector<2, double>& point, const Matrix<2, 2, double>& rotation,
                                  const Vector<2, double>& center) const;
    Vector<2> rotatePoint(const Vector<2>& point, const Matrix<2,2>& rotation) const;
    bool isInsideExcludingRectangle(double angleFrom, double angleTo, double x, double y) const;
    bool checkTangents(double x, double y, double cx, double cy, double r, double angle1, double angle2) const;
    void setAngle(double angle) override;
    bool isIntersection(double p0_x, double p0_y, double p1_x, double p1_y, double p2_x, double p2_y, double p3_x,
                        double p3_y) const;

protected:
    void setPosition(const Vector<2> &position) override;

public:
    static void initClass(const std::string &args);

    Rectangle();
    ~Rectangle() override = default;

    bool overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const override;
    double getVolume(unsigned short dim) const override;
    bool pointInside(BoundaryConditions<2> *bc, const Vector<2> &da) const override;
    bool pointInside(BoundaryConditions<2> *bc, const Vector<2> &da, double angleFrom, double angleTo) const override;

    std::string toWolfram() const override;
    std::string toPovray() const override;
    std::string toString() const override;

    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    Shape<2, 1> *clone() const override;

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;
};


#endif //RSA3D_RECTANGLE_H
