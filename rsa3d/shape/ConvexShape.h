//
// Created by PKua on 12.06.18.
//

#ifndef RSA3D_CONVEXSHAPE_H
#define RSA3D_CONVEXSHAPE_H

#include "Shape.h"

/**
 * @brief A convex shape, which simplifies implementing Shape::voxelInside to implementing ConvexShape::pointInside.
 *
 * It is supported by the fact that it sufficient to check weather all vertices of a voxel lie in Shape's exclusion
 * zone to assure that whole voxel lies in it.
 *
 * @tparam SPATIAL_DIMENSION how many dimensions has the space in which a shape lives
 * @tparam ANGULAR_DIMENSION how many orientational degrees of freedom a shape has. However, this parameter should be
 * non-zero only if ConvexShape::pointInside is capable of dealing with angle-dependent exclusion zones, ie. when
 * generating saturated RSA packings is supported
 */
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
class ConvexShape : public Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> {
protected:
    using EarlyRejectionResult = typename Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::EarlyRejectionResult;

public:
    virtual bool voxelInside(BoundaryConditions<SPATIAL_DIMENSION> *bc, const Vector<SPATIAL_DIMENSION> &voxelPosition,
                     const Orientation<ANGULAR_DIMENSION> &orientation, double spatialSize,
                     double angularSize) const;

    /**
     * @brief Checks if a virtual particle of the same size is within intersection of all excluded volumes for all its
     * orientations denoted by @a orientation and @a orientationRange parameters.
     *
     * @param bc boundary conditions to take into account
     * @param position position of a virtual particle of the same size to check
     * @param orientation array of beginnings of angle intervals
     * @param orientationRange array of lengths of angle intervals
     * @return false if point is outside, true otherwise
     */
    virtual bool pointInside(BoundaryConditions<SPATIAL_DIMENSION> *bc, const Vector<SPATIAL_DIMENSION> &position,
                             const Orientation<ANGULAR_DIMENSION> &orientation,
                             double orientationRange) const = 0;

    /**
     * @brief Checks if a virtual particle of the same size is within excluded volume for any orientation of it.
     *
     * @param bc boundary conditions to take into account
     * @param da position of a virtual particle of the same size to check
     * @return false if point is outside, true otherwise
     */
    virtual bool pointInside(BoundaryConditions<SPATIAL_DIMENSION> *bc,
                             const Vector<SPATIAL_DIMENSION> &position) const;
};

using RSAConvexShape = ConvexShape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>;

#include "ConvexShape.tpp"

#endif //RSA3D_CONVEXSHAPE_H
