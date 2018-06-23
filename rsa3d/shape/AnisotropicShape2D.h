//
// Created by PKua on 07.02.18.
//

#ifndef RSA3D_ANISOTROPICSHAPE2D_H
#define RSA3D_ANISOTROPICSHAPE2D_H


#include "../Matrix.h"
#include "ConvexShape.h"

/**
 * @brief 2D Shape with anisotropy - its position is described by 2 numbers and orientation by 1 number.
 *
 * This class specializes a bit more general Shape class and acts as a base for all anisotropic 2D shapes. Derived
 * classes can now operate on angles, not arrays of angles and automatic normalizataion is introduced (see setAngle()).
 *
 * Note, that Shape::pointInside(BoundaryConditions*,double*,double*,double) and Shape::setOrientation(double*) now
 * delegate to AnisotropicShape2D::pointInside(BoundaryConditions*,double*,double,double) and
 * AnisotripicShape2D::setAngle(double). See those methods' documentation for more information.
 */
class AnisotropicShape2D : public ConvexShape<2, 1> {

protected:

    /**
     * @brief Returns rotation matrix for - getAngle() angle.
     * @return rotation matrix for - getAngle() angle
     */
    Matrix<2, 2> getAntiRotationMatrix() const;

    /**
     * @brief Returns rotation matrix for getAngle() angle.
     * @return rotation matrix for getAngle() angle
     */
    Matrix<2, 2> getRotationMatrix() const;

    /**
     * @brief Add/subrtacts integer multiple of @a interval to @a angleFrom and @a angleTo so that @a angleFrom lies in
     * (0, @a interval) range.
     * @param angleFrom angle range beginning
     * @param angleTo angle range ending
     * @param interval interval to fit angleFrom into
     */
	void normalizeAngleRange(double *angleFrom, double *angleTo, double interval) const;

    /**
     * @brief Add/subrtacts integer multiple of @a interval to @a angle so that it lies in (0, @a interval) and returns
     * result.
     * @param angle angle to normalize
     * @param interval interval to fit angle into
     * @return normalized angle
     */
	double normalizeAngle(double angle, double interval) const;

    /**
     * @brief Delegated to AnisotropicShape2D::setAngle(double). <strong>Override that instead.</strong>
     */
    void setOrientation(const std::array<double, 1> &orientation) final;

    /**
     * @brief Sets new orientation of a shape. It replaces Shape::setOrientation(double*).
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

    /**
     * @brief Delegated to AnisotropicShape2D::pointInside(BoundaryConditions*,double*,double,double). <strong>Implement
     * that instead.</strong>
     */
    bool pointInside(BoundaryConditions<2> *bc, const Vector<2> &position, const std::array<double, 1> &orientation,
                    double orientationRange) const final;

    /**
     * @brief Returns the orientation of a shape.
     * @return the orientation of a shape
     */
    double getAngle() const;

    /**
     * @brief Checks if a given point @a da is within intersection of all excluded zone for all orientations in
     * (@a angleFrom, @a angleTo) range. It replaces Shape::pointInside(BOundaryConditions*,double*,double*,double).
     * @param bc boundary conditions to take into account
     * @param da position of a point to check
     * @param angleFrom array of beginnings of angle intervals
     * @param angleTo array of lengths of angle intervals
     * @return 0 if point is outside, nonzero number otherwise
     */
    virtual bool pointInside(BoundaryConditions<2> *bc, const Vector<2> &da, double angleFrom, double angleTo) const = 0;

};


#endif //RSA3D_ANISOTROPICSHAPE2D_H
