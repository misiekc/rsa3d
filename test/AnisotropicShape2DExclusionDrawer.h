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

    Polygon zone_for_angle(AnisotropicShape2D &shape, double angle, std::size_t resolution);
    Polygon zone_for_two_angles(AnisotropicShape2D &shape, double angle1, double angle2, std::size_t resolution);
    Polygon zone_for_angle_range(AnisotropicShape2D &shape, double angleFrom, double angleTo,
                                 std::size_t resolution);
    void polygon_to_wolfram(const Polygon &polygon, std::ostream &stream);
    void print_notebook(AnisotropicShape2D &shape, const Polygon &fromZone, const Polygon &toZone,
                        const Polygon &fromAndToZone, const Polygon &rangeZone, std::ostream &stream);
    int main(int argc, char **argv);
}

#endif //RSA3D_ANISOTROPICSHAPE2DEXCLUSIONDRAWER_H
