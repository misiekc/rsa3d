//--------------------------------------------------------------------------------------------
// Test of cuboid intersection algorithm. It uses standard triangle-triangle intersection
// algorithm to check results given by Cuboid::overlap
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _INTTEST_H
    #define _INTTEST_H
    
#include <ostream> 
#include <vector>
#include "../rsa3d/Parameters.h"
#include "../rsa3d/shapes/OverlapStrategy.h"
#include "utility/ShapePairFactory.h"

/**
 * @brief A test of Shape::overlap method.
 *
 * Usage:
 * <blockquote>./rsa_test shape_ovtest [particle] [attibutes] [ball_radius] [max_tries]</blockquote>
 *
 * The test is using BallFactory for creating shape pairs. It can be performed provided a shape configured in input.txt
 * is a subclass of OverlapStrategyShape and supports at least two overlap strategies.
 *
 * Each second and following strategies will be compared with the first one (the order from
 * OverlapStrategyShape::getStrategies) and each pair with two different overlap results will be dump to ovtest_dump.nb
 * file.
 */
namespace shape_ovtest
{
    /**
     * @brief Main function for the test of Shape::overlap method.
     *
     * Usage:
     * <blockquote>./rsa_test shape_ovtest [particle] [attibutes] [ball_radius] [max_tries]</blockquote>
     */
    int main(int argc, char **argv);
}

#endif  // _INTTEST_H
