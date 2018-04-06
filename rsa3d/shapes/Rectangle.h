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

    // from: http://geomalgorithms.com/a03-_inclusion.html

    /* isLeft(): tests if a point is Left|On|Right of an infinite line.
     * Input:  three points P0, P1, and P2
     * Return:
     * >0 for P2 left of the line through P0 and P1
     * =0 for P2  on the line
     * <0 for P2  right of the line
     */
    double isLeft(double p0x, double p0y, double p1x, double p1y, double p2x, double p2y) const;

    /*
     * wn_PnPoly(): winding number test for a point in a polygon
     * Input:   x,y = a point,
     * xs, ys = vertex points of a polygon xs[n+1] with xs[n]=xs[0]
     * Return:  wn = the winding number (=0 only when P is outside)
     */
    int wn_PnPoly(double x, double y, const double *xs, const double *ys, int n) const;

protected:
    void setPosition(const double *position) override;

public:
    static void initClass(const std::string &args);
    static Shape<2, 1> * create(RND *rnd);

    Rectangle();
    Rectangle(const Rectangle &other);
    Rectangle & operator = (const Rectangle & el);
    ~Rectangle() override = default;

    int overlap(BoundaryConditions *bc, Shape<2, 1> *s) const override;

    double getVolume() const override;

    double getVoxelAngularSize() const override;

    int pointInside(BoundaryConditions *bc, double* da) const override;
    int pointInside(BoundaryConditions *bc, double* da, double angleFrom, double angleTo) const override;

    double getNeighbourListCellSize() const override;
    double getVoxelSize() const override;


    std::string toWolfram() const override;
    std::string toPovray() const override;
    std::string toString() const override;

    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
};


#endif //RSA3D_RECTANGLE_H
