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
    static double neighbourListCellSize;
    static double voxelSize;
    static double longer;
    static double shorter;

    // b is shorter
    double a, b, halfA, halfB;

    double getXOnCircle(double cx, double cy, double r, double angle);

    double getYOnCircle(double cx, double cy, double r, double angle);

    int isInsideCircle(const Vector<2, double> &point, double cx, double cy, double r) const;

    int isInsideRect(const Vector<2, double> &point, double cx, double cy, double a, double b);

    Vector<2, double> rotatePoint(const Vector<2, double>& point, const Matrix<2, 2, double>& rotation, const Vector<2, double>& center);

    Vector<2> rotatePoint(const Vector<2>& point, const Matrix<2,2>& rotation);

    void getNewSize(const Matrix<2, 2, double> &rotationFrom, const double *p, const Vector<2, double> &newSizeFrom) const;

    bool checkExcludingRectangle(double angleFrom, double angleTo, double x, double y) const;

    bool checkTangents(double x, double y, double cx, double cy, double r, double angle1, double angle2);

    void setAngle(double angle);

public:
    static void initClass(const std::string &args);
    static Shape<2, 1> * create(RND *rnd);

    Rectangle();
    Rectangle(const Rectangle &other);
    Rectangle & operator = (const Rectangle & el);
    ~Rectangle() override = default;

    void rotate(double *v) override;

    int overlap(BoundaryConditions *bc, Shape<2, 1> *s) override;

    double getVolume() override;

    int pointInside(BoundaryConditions *bc, double* da) override;
    int pointInside(BoundaryConditions *bc, double* da, double angleFrom, double angleTo) override;

    double getNeighbourListCellSize() override;
    double getVoxelSize() override;


    std::string toWolfram() const override;
    std::string toPovray() const override;
    std::string toString() override;

    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
};


#endif //RSA3D_RECTANGLE_H
