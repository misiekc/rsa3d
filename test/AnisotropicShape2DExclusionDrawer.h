//
// Created by PKua on 31.03.18.
//

#ifndef RSA3D_ANISOTROPICSHAPE2DEXCLUSIONDRAWER_H
#define RSA3D_ANISOTROPICSHAPE2DEXCLUSIONDRAWER_H

#include <string>
#include <vector>
#include "../rsa3d/AnisotropicShape2D.h"
#include "../rsa3d/Vector.h"


/**
 * @brief A simple program generating Wolfram Mathematica notebook presenting exclusion zone for AnisotropicShape2D.
 *
 * Usage:
 * <blockquote>
 * rsa as2d_exdrawer (shape) (attr) (shape angle) (angle from) (angle to) (resolution) (output file)
 * </blockquote>
 *
 * The output file contains:
 * <ul>
 * <li> exclusion zones for angle from and angle to
 * <li> intersection on two above zones
 * <li> intersection of all ranges from given angle range (Shape::pointInside)
 * </ul>
 *
 * It can be used after checking exclusion zones correctness using as2d_extest (because it assumes that regions are
 * convex) for drawings in an article or more accurate examination of a shape of those zones.
 */
namespace as2d_exdrawer
{
    /**
     * More descriptive name for the vector of vertices of a polygon.
     */
    using Polygon = std::vector<Vector<2>>;

    /**
     * Returns a as2d_exdrawer::Polygon describing exclusion zone for given @a shape and @a angle, sampling @a
     * resolution number of points and testing Shape::overlap()
     * @param shape shape which exclusion zone to calculate
     * @param angle angle of exclusion zone
     * @param resolution number of vertices of resulting Polygon
     * @return as2d_exdrawer::Polygon describing exclusion zone
     */
    Polygon zone_for_angle(AnisotropicShape2D &shape, double angle, std::size_t resolution);

    /**
     * Same as zone_for_angle(), but resulting as2d_exdrawerPolygon is intersection of exclusion zones for @a angle1 and
     * @a angle2
     */
    Polygon zone_for_two_angles(AnisotropicShape2D &shape, double angle1, double angle2, std::size_t resolution);

    /**
     * Same as zone_for_angle(), but resulting as2d_exdrawer::Polygon is exclusion zone for angle range, sampled using
     * Shape::pointInside
     */
    Polygon zone_for_angle_range(AnisotropicShape2D &shape, double angleFrom, double angleTo,
                                 std::size_t resolution);

    /**
     * Converts given as2d_exdrawer::Polygon to Wolfram format and prints it on stream
     */
    void polygon_to_wolfram(const Polygon &polygon, std::ostream &stream);

    /**
     * Creates a whole Wolfram notebook visualizing exclusion zone for given @a shape and prints it on @a stream
     * @param shape shape, drawn Black, which exclusion zone is visualized
     * @param fromZone exclusion zone for fromAngle, drawn Red
     * @param toZone exclusion zone for toAngle, drawn Yellow
     * @param fromAndToZone intersection of from- and toZone, drawn Purple; it sticks from behind rangeZone
     * @param rangeZone exclusion zone for angle range, drawn Orange
     * @param stream stream to print on
     */
    void print_notebook(AnisotropicShape2D &shape, const Polygon &fromZone, const Polygon &toZone,
                        const Polygon &fromAndToZone, const Polygon &rangeZone, std::ostream &stream);

    /**
     * Main function for exclusion drawind facility
     *
     * Usage:
     * <blockquote>
     * rsa as2d_exdrawer (shape) (attr) (shape angle) (angle from) (angle to) (resolution) (output file)
     * </blockquote>
     */
    int main(int argc, char **argv);
}

#endif //RSA3D_ANISOTROPICSHAPE2DEXCLUSIONDRAWER_H
