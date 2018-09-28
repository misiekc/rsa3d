//--------------------------------------------------------------------------------------------
// Test of Cuboid::pointInside method
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CUBOID_POINT_INSDE_TEST_H
    #define _CUBOID_POINT_INSDE_TEST_H

/**
 * @brief A test of Shape::pointInside method.
 *
 * Usage:
 * <blockquote>./rsa_test shape_pitest [particle] [attibutes] [ball_radius] [max_tries]</blockquote>
 *
 * The test is using BallFactory for creating shape pairs. For each pair it checks results given by Shape::overlap and
 * Shape::pointInside. If Shape::pointInside gives true, so should Shape::overlap do - otherwise a @a conflict occurs.
 */
namespace shape_pitest
{
    /**
     * @brief Main function for the test of Shape::pointInside method.
     *
     * Usage:
     * <blockquote>./rsa_test shape_pitest [particle] [attibutes] [ball_radius] [max_tries]</blockquote>
     */
    int main(int argc, char **argv);
}

#endif // _CUBOID_POINT_INSDE_TEST_H
