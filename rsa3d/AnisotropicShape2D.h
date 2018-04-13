//
// Created by PKua on 07.02.18.
//

#ifndef RSA3D_ANISOTROPICSHAPE2D_H
#define RSA3D_ANISOTROPICSHAPE2D_H


#include "Shape.h"
#include "Matrix.h"

class AnisotropicShape2D : public Shape<2, 1> {

protected:

    /**
     * Returns rotation matrix for -getAngle() angle.
     * @return rotation matrix for -getAngle() angle
     */
    Matrix<2, 2> getAntiRotationMatrix() const;

    /**
     * Returns rotation matrix for getAngle() angle.
     * @return rotation matrix for getAngle() angle
     */
    Matrix<2, 2> getRotationMatrix() const;

    /**
     * Add/subrtacts integer multiple of an interval to angleFrom and angleTo so that angleFrom lies in (0, interval)
     * range.
     * @param angleFrom angle range beginning
     * @param angleTo angle range ending
     * @param interval interval to fit angleFrom into
     */
	void normalizeAngleRange(double *angleFrom, double *angleTo, double interval) const;

    /**
     * Add/subrtacts integer multiple of an interval to angle so that it lies in (0, interval) and returns result.
     * @param angle angle to normalize
     * @param interval interval to fit angle into
     * @return normalized angle
     */
	double normalizeAngle(double angle, double interval) const;

    // setOrientation(const double*) is delegated to setAngle(double)
    void setOrientation(const double *orientation) final;

    /**
     * Sets new orientation of a shape. It replaces Shape::setOrientation(double*).
     *
     * The angle is normalized using normalizeAngle(double, double) with getVoxelAngularSize() as an interval.
     *
     * Derived shape classes have to override this method when they want to keep track of shape's orientation. One would
     * then typically write:
     * \code
     * void Derived::setAngle(double angle) {
     *     AnisotropicShape2D::setAngle(angle);
     *     // invoke getAngle and for example compute vertices
     * }
     * \endcode
     * @param angle new orientation of a shape
     */
    virtual void setAngle(double angle);

public:

    // pointInside(double* orientation, double orientationRange) -> pointInside(double angleFrom, double angleTo)
    int pointInside(BoundaryConditions *bc, double* position, double *orientation, double orientationRange) const final;

    /**
     * Returns the orientation of a shape.
     * @return the orientation of a shape
     */
    double getAngle() const;

    /**
     * Checks if a given point is within intersection of all excluded zone for all orientations in (angleFrom, angleTo)
     * range.
     * @param bc boundary conditions to take into account
     * @param da position of a point to check
     * @param angleFrom array of beginnings of angle intervals
     * @param angleTo array of lengths of angle intervals
     * @return 0 if point is outside, nonzero number otherwise
     */
    virtual int pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) const = 0;
};


#endif //RSA3D_ANISOTROPICSHAPE2D_H
