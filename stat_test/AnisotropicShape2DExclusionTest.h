//
// Created by PKua on 08.02.18.
//

#ifndef RSA3D_SPHEROCYLINDER2DINTTEST_H
#define RSA3D_SPHEROCYLINDER2DINTTEST_H

#include "../rsa3d/geometry/Vector.h"
#include "../rsa3d/shape/AnisotropicShape2D.h"
#include "utils/ShapePairFactory.h"
#include <memory>

/**
 * @brief A simple program generating Wolfram Mathematica notebook for quick testing of exclusion zone for
 * AnisotropicShape2D.
 *
 * Usage:
 * <blockquote>
 * ./rsa_test as2d_extest <input> <output = out.nb>
 * </blockquote>
 *
 * It fixes the first shape of a given position and orientation and checks the results of the one of:
 * <ul>
 * <li> Shape::pointInside with different random points in the area of the first shape
 * <li> Shape::overlap with the second shape of fixed orientation but different random positions in the area of the
 * first shape
 * </ul>
 * Each checked point is drawn in resulting Mathematica notebook in appropriate color as described later. Shapes are
 * also printed if a shape class provides wolfram drawing facility. These points form exclusion zones which can be then
 * verified by eye. Note that, unlike as2d_exdrawer, these zones can be concave or totally messed up, so it is the test
 * better suited for debugging.
 *
 * The input file is configuration file with lines in a form of key=value containing keys:
 * <ul>
 * <li> <strong>particle</strong> - ShapeFactory name of a particles to test - WORKS ONLY FOR Shape<2, 1>
 * <li> <strong>attr</strong> - ShapeFactory attributes of a particle
 * <li> <strong>mode</strong> - @a overlap or @a point_inside, depending on what test to perform
 * <li> <strong>first_x</strong> - x position of the first shape
 * <li> <strong>first_y</strong> - y position of the first shape
 * <li> <strong>first_angle</strong> - angle of the first shape IN DEGREES
 * <li> <strong>box_width</strong> - width of a box from which random positions will be chosen (the first shape lies in
 * the middle)
 * <li> <strong>box_height</strong> - height of a box
 * <li> <strong>points</strong> - number of points to generate and check
 * <li> <strong>append</strong> - @a true of @a false. If true, the resulting notebook will be appended to output file
 * </ul>
 *
 * Moreover, for given modes there are additional fields.
 * <ul>
 *     <li>
 *     For @a overlap:
 *     <ul>
 *         <li> <strong>second_angle</strong> - the orientation of the second shape IN DEGREES
 *     </ul>
 *     </li>
 *     <li>
 *     For @a point_inside:
 *     <ul>
 *         <li> <strong>from_angle</strong>, <strong>to_angle</strong> - angle range for ConvexShape::point_inside
 *     </ul>
 *     </li>
 * </ul>
 *
 * And finally, colors of points:
 * <ul>
 * <li> <strong>Black</strong> - color of the first shape
 * <li> <strong>Blue</strong> - (for @a overlap) color of the second shape
 * <li> <strong>Orange</strong> - color of points inside exclusion zone (the type of a zone depends on test mode)
 * <li> <strong>Green</strong> - color of points outside all "interesting" areas
 * <li> <strong>Purple</strong> - (for @a point_inside) color of points in both Shape::overlap exclusion zones for
 * endpoints of ConvexShape::pointInside angle range but eventually not in a ConvexShape::pointInside exclusion zone
 * <li> <strong>Red</strong> - (for @a point_inside) color of points in a Shape::overlap exclusion zone for the first
 * angle of ConvexShape::pointInside but not in the second exclusion zone
 * <li> <strong>Yellow</strong> - (for @a point_inside) color of points in a Shape::overlap exclusion zone for the
 * second angle of ConvexShape::pointInside but not in the first exclusion zone
 * </ul>
 */
namespace as2d_extest
{
    bool point_inside(const Shape<2, 1> &shape, const Vector<2> &point, double angleFrom, double angleTo);

    /**
     * @brief Main function for exclusion testing facility.
     *
     * Usage:
     * <blockquote>
     * ./rsa_test as2d_extest <input> <output = out.nb>
     * </blockquote>
     */
    int main(int argc, char **argv);
}

#endif //RSA3D_SPHEROCYLINDER2DINTTEST_H
