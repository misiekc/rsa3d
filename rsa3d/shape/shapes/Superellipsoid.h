/*
* Created by Michal Ciesla on 30.06.2022.
*
* Author: Michal Ciesla
* Sources:
*      	[2] "Jammed Packings of Hard Particles" Aleksandar Donev, 2006
*
*/

#ifndef RSA3D_SUPERELLIPSOID_H
#define RSA3D_SUPERELLIPSOID_H

#include "../ConvexShape.h"
#include "../OrderCalculable.h"

class Superellipsoid : public ConvexShape<3, 0>, public OrderCalculable {
private:
    static void normalizeVolume();
    static double a, b, c, p;

    Matrix<3, 3> orientation;		// Cartesian rotation matrix 3D
    Matrix<3, 3> orientationTr;                // orientation matrix transposed
    Matrix<3, 3> orientationInv;                // orientation matrix inversed

    Vector<3> r_relative(const Vector<3> &r) const;
    double inner_shape_function(const Vector<3> &r) const;
    double shape_function(const Vector<3> &r) const;
    Vector<3> inner_shape_function_gradient(const Vector<3> &r) const;
    Vector<3> shape_function_gradient(const Vector<3> r) const;
    Matrix<3, 3> inner_shape_function_hesian(const Vector<3> r) const;
    Matrix<3, 3> shape_function_hesian(const Vector<3> r) const;

    Vector<3> gradient(const Vector<3> r) const;
    Vector<3> nabla(const Vector<3> r) const;
    Matrix<3,3> hesian(const Vector<3> r) const;
    Matrix<3,3> nabla_2(const Vector<3> r) const;

public:
    static void initClass(const std::string &attr);

    explicit Superellipsoid(const Matrix<3, 3> &orientation);



    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const override;
    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     double orientationRange) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    Shape<3, 0> *clone() const override;

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;

    std::string toPovray() const;
    std::string toWolfram() const;
//    Matrix<3, 3> getEllipsoidMatrix() const;

};

#endif //RSA3D_ELLIPSOID_H
