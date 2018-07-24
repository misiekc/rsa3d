/*
* Created by PKua on 17.07.18.
*
* Author: Wojciech Nowak
* Sources:
*       [1] "Statistical Mechanics of Hard Ellipsoids. I. Overlap Algorithm and the Contact Function", John W. Perram, M. S. Wertheim, February 9, 1984
*      	[2] "Jammed Packings of Hard Particles" Aleksandar Donev, 2006
*
*/

#ifndef RSA3D_ELLIPSOID_H
#define RSA3D_ELLIPSOID_H

#include "../ConvexShape.h"

const double INITIAL_GUESS = 0.5;	// suggested by A.Donev [2] initial guess evaluation is not applicable in our case as we have the same ellipsoids

class Ellipsoid : public ConvexShape<3, 0> {
private:
    static void normalizeVolume();
    static double a, b, c;

    Matrix<3, 3> orientation;		// Cartesian rotation matrix 3D
    Matrix<3, 3> X;                 // Matrix with axis length on diagonal (2.5)[1]
    Matrix<3, 3> M;                 // Ellipsoid specific matrix (2.5)[1]

    inline Matrix<3, 3> matrixC(const Matrix<3, 3> &A, const Matrix<3, 3> &B, double lambda) const;
    inline double quadraticForm(const Matrix<3, 3> &M, const Vector<3> &v) const;
    inline Vector<3> vectorRAB(const Vector<3> & rA, const Vector<3> & rB) const;
    inline double overlapFunction(const Matrix<3, 3> &A, const Matrix<3, 3> &B, const Vector<3> &rAB, double lambda) const;
    inline double firstDerivative(const Matrix<3, 3> &A, const Matrix<3, 3> &B, const Vector<3> &rAB, double lambda) const;
    inline double secondDerivative(const Matrix<3, 3> &A, const Matrix<3, 3> &B, const Vector<3> &rAB, double lambda) const;
    double getNewLambda(double lambda, double derivativeQuotient) const;

public:
    static void initClass(const std::string &attr);

    explicit Ellipsoid(const Matrix<3, 3> &orientation) : orientation(orientation) {}

    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const override;
    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     double orientationRange) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    Shape<3, 0> *clone() const override;

    std::string toPovray() const;
    Matrix<3, 3> getEllipsoidMatrix() const;
};

#endif //RSA3D_ELLIPSOID_H