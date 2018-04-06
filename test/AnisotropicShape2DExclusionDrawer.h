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
 * A simple program generating Wolfram Mathematica notebook presenting exclusion zone for AnisotropicShape2D.
 *
 * Usage: rsa as2d_exdrawer (shape) (attr) (shape angle) (angle from) (angle to) (resolution) (output file)
 *
 * The output file contains:
 * <ul>
 * <li> exclusion zones for angle from and angle to
 * <li> intersection on two above zones
 * <li> intersection of all ranges from given angle range (pointInside)
 * </ul>
 */
namespace as2d_exdrawer
{
    using Polygon = std::vector<Vector<2>>;

    /**
     * Returns a Polygon describing exclusion zone for given shape and angle, sampling resolution points and testing
     * overlap function
     * @param shape shape which exclusion zone to calculate
     * @param angle angle of exclusion zone
     * @param resolution number of vertices of resulting Polygon
     * @return Polygon describing exclusion zone
     */
    Polygon zone_for_angle(AnisotropicShape2D &shape, double angle, std::size_t resolution);

    /**
     * Same as zone_for_angle, but resulting Polygon is intersection of exclusion zones for angle1 and angle2
     */
    Polygon zone_for_two_angles(AnisotropicShape2D &shape, double angle1, double angle2, std::size_t resolution);

    /**
     * Same as zone_for_angle, but resulting Polygon is exclusion zone for angle range, sampled using pointInside
     */
    Polygon zone_for_angle_range(AnisotropicShape2D &shape, double angleFrom, double angleTo,
                                 std::size_t resolution);

    /**
     * Converts given polygon to Wolfram format and prints it on stream
     */
    void polygon_to_wolfram(const Polygon &polygon, std::ostream &stream);

    /**
     * Creates a whole Wolfram notebook visualizing exclusion zone for given shape and prints it on stream
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
     * Usage: rsa as2d_exdrawer (shape) (attr) (shape angle) (angle from) (angle to) (resolution) (output file)
     */
    int main(int argc, char **argv);
}

#endif //RSA3D_ANISOTROPICSHAPE2DEXCLUSIONDRAWER_H
