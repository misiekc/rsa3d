//
// Created by PKua on 12.08.18.
//

#ifndef RSA3D_SHAPEBCTEST_H
#define RSA3D_SHAPEBCTEST_H

#endif //RSA3D_SHAPEBCTEST_H

/**
 * @brief A test of correct usage of BoundaryConditions in Shape::overlap and ConvexShape::pointInside.
 *
 * Usage:
 * <blockquote>./rsa_test shape_bctest [particle] [attibutes] [ball_radius] [max_tries]</blockquote>
 *
 * The test generates random pairs uniformly from a ball and random translations from the same point distribution. Then,
 * it compares results given by Shape::overlap and ConvexShape::pointInside passing BoundaryConditions giving
 * beforementioned translation with those methods evaluated for a pair with translation already applied. If results
 * differ, a conflict occurs and is reported.
 */
namespace shape_bctest
{
    /**
     * @brief A main function for shape_bctest.
     *
     * Usage:
     * <blockquote>./rsa_test shape_bctest [particle] [attibutes] [ball_radius] [max_tries]</blockquote>
     */
    int main(int argc, char **argv);
}